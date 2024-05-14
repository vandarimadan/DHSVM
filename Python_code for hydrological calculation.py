# Title: Creating stream Network and Stream Map for DHSVM Model
# ORG:Pacific Northwest National Laboratory
# AUTHOR:Zhuoran Duan
# ORIG-DATE:    Apr-2017
# Updated by: Madan Raj Bhandari
#Organization: Memorial University of Newfoundland, Canada
#Position: Graduate student, Department of Biology
# COMMENTS: This python script is modified based on original 
#-AML scripts createstreamnetwork.aml and updated python script developed by Zhuoran Duan for running the script to particular watershed as part of DHSVM
# Last Change: 2024-05-09 

#Import modules for creating stream network and stream Map

import arcpy
from arcpy import env
from arcpy.sa import *
import arcgisscripting
import csv 
import math
import sys
import os
import numpy as np

#Defining folder or workspace location

env.workspace = "D:\\calibration\\Python_McKinley\\"   
path = "D:\\calibration\\Python_McKinley\\"   

#initial inputs setup and defining minimum contributing area and min and max soil depth

dem = "DEM"                         #input the name of dem on "....." (elevation raster)
wshed = "Watershed"                 #input the name of watershed boundary raster file or mask raster delinating watershed boundary
soildepth = "soil_depth"            # name for the soil depth file that will be generated based on min and max soil depth
mindepth =2.5                       # minimum soil depth (meter)
maxdepth = 5.0                      # maximun soil depth (meter)
streamfile = "streamfile"           #name of stream file that will be generated
con_area = 10000000                 #minimum contributing area in meter square (larger the contributing area lesser that stream segment)
key = 'Mask'                        # Mask or point; input mask if you have watershed boundary raster or input point if watershed will be generated from point feature

#Ensure whether your GIS platform has spatial analysis extension or not

#Verifying the input DEM and watershed raster
if not arcpy.Exists(dem):
    sys.exit("Digital Elevation Model (DEM) input not valid, exit program")

if not arcpy.Exists(wshed):
    sys.exit("watershed/point data input not valid, exit program")

#setting the cell size and extension for the analysis; cell size (i.e. resolution) and extension (i.e. number of rows and columns should be same) 
env.outputCoordinateSystem = dem
env.extent = dem
env.cellSize = dem
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

#creating a watershed mask raster of full extension based on dem extension
if arcpy.Exists("watmask"):
    arcpy.Delete_management("watmask")
MaskExt=CreateConstantRaster (1, "INTEGER", env.cellSize, env.extent)
MaskExt.save("watmask")

#generating flow direction raster file from DEM
if arcpy.Exists("flow_dir"):
    arcpy.Delete_management("flow_dir")

print('Flow direction')
flowdir = FlowDirection(dem, "", "")

# Generating watershed if the key is "point" 
if key=='Mask':
    print('Mask/watershed boundary file provided')
    wshd = wshed
    env.mask = wshd
    flowdir = FlowDirection(dem, "", "")
elif key=='MOUTH':
    print('Mask/watershed boundary file not provided, generating watershed')
    flowdir = FlowDirection(dem, "", "")
    wshd = Watershed(flowdir , wshed, "")
    env.mask = wshd
else:
    print('Input key not valid')

flowdir.save("flow_dir")

# Computing flow accumulation
if arcpy.Exists("flowacc"):
    arcpy.Delete_management("flowacc")
env.mask = wshd    
print('Generating flow accumulation')    
flowacc = FlowAccumulation(flowdir, "", "INTEGER")
flowacc.save("flowacc")

detlaxResult = arcpy.GetRasterProperties_management(dem,"CELLSIZEX")
deltax = detlaxResult.getOutput(0)
sourcepix = con_area / (float(deltax) * float(deltax))          #contributing area (con_area) defined above is used here to generate stream network

temp = Con(flowacc > sourcepix, 1)  

rivg= temp * wshd

#Creating Geodata base to store data
if arcpy.Exists("/output.gdb"):
    print('geodatabase already exist')
else:
    print('geodatabase does NOT exist, create geodatabse')
    streamgdb = "output.gdb"
    arcpy.CreateFileGDB_management(path, streamgdb)

streamlink=StreamLink (rivg, flowdir)
if arcpy.Exists(path+"/output.gdb/"+streamfile):
    arcpy.Delete_management(path+"/output.gdb/"+streamfile)
    print('deleted previous stream file and create a new one')

streamnetwork=path+"/output.gdb/"+streamfile
print('creating stream shapefile')
StreamToFeature(streamlink, flowdir, streamnetwork, "NO_SIMPLIFY")    

if arcpy.Exists(rivg):
    arcpy.Delete_management(rivg)

#Calculating contributing area for each cell 
arcpy.PolylineToRaster_conversion (streamnetwork, "arcid", 'rast', "MAXIMUM_LENGTH","NONE", env.cellSize)

local= Watershed(flowdir, 'rast', "VALUE")          #creating a sub-catchment based on contributing area

if arcpy.Exists("local"):
    arcpy.Delete_management("local")

local.save("local")
arcpy.AddField_management (streamnetwork, "local", "LONG")
arcpy.JoinField_management(streamnetwork,"arcid",local,"VALUE","COUNT")  #sometime count value of 'local' raster may not be joined and error may occur (value can be enter manually).

fields=['COUNT','local']

with arcpy.da.UpdateCursor(streamnetwork, fields) as cursor:
    for row in cursor:
        if row[0] is None:
            row[1] = 0
        else:
            row[1]=row[0]
        cursor.updateRow(row)
del row

arcpy.DeleteField_management(streamnetwork, "COUNT")

arcpy.Delete_management('rast') #deleting Unwanted file

# creating node point for each stream segment (start node and end node), finding elevations and calculating slope

if arcpy.Exists(path+"output.gdb/nodestart"):
    arcpy.Delete_management(path+"/output.gdb/nodestart")
    print('start node file already exists, deleting and creating new one')

nodestart=path+"output.gdb/nodestart"
arcpy.FeatureVerticesToPoints_management(streamnetwork, nodestart, "START")

if arcpy.Exists(path+"output.gdb/nodeend"):
    arcpy.Delete_management(path+"/output.gdb/nodeend")
    print('end node file already exists, delete and create new')

nodeend=path+"output.gdb/nodeend"
arcpy.FeatureVerticesToPoints_management(streamnetwork, nodeend, "END")

#calculating node start elevation point and nodeend elevation point to calculate slope

env.mask="watmask"
tmpacc=Con(IsNull("flowacc")==1,int(flowacc.maximum),"flowacc") #temporary accumulation value

ExtractMultiValuesToPoints(nodestart, [[tmpacc, "MAXGRID"]])

arcpy.AddField_management (streamnetwork, "downarc", "LONG")

ExtractMultiValuesToPoints(nodestart, [[dem, "SELEV"]], "NONE")
demras=Raster(dem)
minimum_value = demras.minimum
#save demras as a file if error occur
tmpdem = Con(IsNull(demras),minimum_value, dem)
ExtractMultiValuesToPoints(nodeend, [[tmpdem, "EELEV"]], "NONE")
env.mask = wshd  

#adding columns for hydrological calculation

arcpy.JoinField_management(streamnetwork,"arcid",nodestart,"arcid","SELEV")
arcpy.JoinField_management (streamnetwork, "arcid", nodeend, "arcid", "EELEV")
arcpy.JoinField_management (streamnetwork, "arcid", nodestart, "arcid", "MAXGRID")

arcpy.AddField_management (streamnetwork, "uparc", "LONG")
arcpy.AddField_management (streamnetwork, "dz", "FLOAT", 12, 3)
arcpy.AddField_management (streamnetwork, "slope", "FLOAT", 12, 5)
arcpy.AddField_management (streamnetwork, "meanmsq",  "FLOAT")
arcpy.AddField_management (streamnetwork, "segorder", "LONG")
arcpy.AddField_management (streamnetwork, "chanclass", "SHORT", 8, "")
arcpy.AddField_management (streamnetwork, "hyddepth", "FLOAT", 8, 2)
arcpy.AddField_management (streamnetwork, "hydwidth", "FLOAT", 8, 2)
arcpy.AddField_management (streamnetwork, "effwidth", "FLOAT", 8, 2)
arcpy.AddField_management (streamnetwork, "effdepth", "FLOAT", 8, 2)

print('Calculating Slope for each stream segment')
arcpy.CalculateField_management (streamnetwork, "dz", "abs(!SELEV! - !EELEV!)", "PYTHON3", "")

expression = "clacSlope(float(!dz!),float(!Shape_Length!))"   #slope calculation
codeblock = """
def clacSlope(dz, length):
    if (dz / length) > 0.00001:
        return dz / length
    else:
        return 0.00001
"""


arcpy.CalculateField_management (streamnetwork, "slope", expression, "PYTHON3", codeblock)

arcpy.CalculateField_management (streamnetwork, "segorder", "1", "PYTHON3")
arcpy.CalculateField_management (streamnetwork, "uparc", "0", "PYTHON3")
arcpy.CalculateField_management (streamnetwork, "downarc", "-1", "PYTHON3")
arcpy.CalculateField_management (streamnetwork, "meanmsq", "0.0", "PYTHON3")

# Calculating Segment Order  and mean meter square for each sub-basin    
#-------------------------------------------------------------------#
arr=arcpy.da.TableToNumPyArray(streamnetwork,('from_node','to_node','segorder','local','MAXGRID','meanmsq','uparc','downarc','arcid'))

print('Looking for downstream arc')
# Calculate downstream arc
for jj, ii in enumerate(arr['to_node']):
    arr2=[]
    for i, j in enumerate(arr['from_node']):
        if j == ii:
            arr['downarc'][jj]=arr['arcid'][i]


print('Looking for upstream arc')
# Calculate upstream arc 
for jj, ii in enumerate(arr['from_node']):
    arr2=[]
    for i, j in enumerate(arr['to_node']):
        if j == ii:
            arr2=np.append(arr2,i)
  

#input env.cellsize value (in my case its 30.0)here 
if not len(arr2):
        arr['uparc'][jj]=-1
        arr['segorder'][jj]=1
        arr['meanmsq'][jj]=arr['local'][jj]/ 2 * float(30.0) * float(30.0)
else:
        arr3=arr2.astype(int)
        loc=np.argmax(arr['MAXGRID'][arr3]+arr['local'][arr3])
        arr['uparc'][jj]=arr['arcid'][arr3[loc]]
        arr['meanmsq'][jj]=(arr['MAXGRID'][jj]+arr['local'][jj]/ 2) * float(30.0) * float(30.0)

# Calculating segorder for each stream segment
order=1
a=99
print('Calculating segment order')
while a > 0:
    a=0
    for jj, ii in enumerate(arr['segorder']):
        if ii==order:
            a+=1
            for i, j in enumerate(arr['arcid']):
                if j == arr['downarc'][jj]:
                    arr['segorder'][i]=max(order+1,arr['segorder'][i])            
    order+=1
   
               
# inserting the array to streamnetwork table
arcpy.da.ExtendTable(streamnetwork,  "arcid", arr, "arcid", append_only=False)

# if the above code doesnot work to calculate meanmsq; use these code
fields = ['local', 'meanmsq']

# Using UpdateCursor to iterate and update the 'meanmsq' field
with arcpy.da.UpdateCursor(streamnetwork, fields) as cursor:
    for row in cursor:
        # row[0] is the value in 'local', row[1] is the value in 'meanmsq'
        # Perform the calculation: local / 2 * 30.0 * 30.0
        row[1] = (row[0] / 2) * (30.0 ** 2)  # 30 squared is 900, so you multiply local/2 by 900
        cursor.updateRow(row)

print("Update complete.")


# stream segment Hydraulic properties and stream class       
import sys
import arcpy
from arcpy import env
from arcpy.sa import *
import arcgisscripting

 # row[0] - 'slope'
        # row[1] - 'meanmsq' - mean contributing area to stream segment
        # row[2] - 'chanclass' - stream class number
        # row[3] - 'hydwidth'
        # row[4] - 'hyddepth'
        # row[5] - 'effwitdth'
fields=['slope','meanmsq','chanclass','hydwidth','hyddepth','effwidth']

with arcpy.da.UpdateCursor(streamnetwork, fields) as cursor:
    for row in cursor:
            if (row[0] <= 0.002 and row[1] <= 1000000):
                row[2] = 1
                row[3] = 0.5
                row[4] = 0.03
                row[5] = 0.06
            elif (row[0] <= 0.002 and (row[1] > 1000000 and row[1] <= 10000000)):
                row[2] = 2
                row[3] = 1.0
                row[4] = 0.03
                row[5] = 0.09
            elif (row[0] <= 0.002 and (row[1] > 10000000 and row[1] <= 20000000)):
                row[2] = 3
                row[3] = 2.0
                row[4] = 0.03
                row[5] = 0.12
            elif (row[0] <= 0.002 and (row[1] > 20000000 and row[1] <= 30000000)):
                row[2] = 4
                row[3] = 3.0
                row[4] = 0.03
                row[5] = 0.15
            elif (row[0] <= 0.002 and (row[1] > 30000000 and row[1] <= 40000000)):
                row[2] = 5
                row[3] = 4.0
                row[4] = 0.03
                row[5] = 0.18
            elif (row[0] <= 0.002 and row[1] > 40000000) :
                row[2] = 6
                row[3] = 4.5
                row[4] = 0.03
                row[5] = 0.21
            elif ((row[0] > 0.002 and row[0] <= 0.1) and row[1] <= 1000000):
                row[2] = 7
                row[3] = 0.5
                row[4] = 0.05
                row[5] = 0.1
            elif ((row[0] > 0.002 and row[0] <= 0.1) and (row[1] > 1000000 and row[1] <= 10000000)):
                row[2] = 8
                row[3] = 1.0
                row[4] = 0.05
                row[5] = 0.15
            elif ((row[0] > 0.002 and row[0] <= 0.1) and (row[1] > 10000000 and row[1] <= 20000000)):
                row[2] = 9
                row[3] = 2.0
                row[4] = 0.05
                row[5] = 0.2
            elif ((row[0] > 0.002 and row[0] <= 0.1) and (row[1] > 20000000 and row[1] <= 30000000)):
                row[2] = 10
                row[3] = 3.0
                row[4] = 0.05
                row[5] = 0.25
            elif ((row[0] > 0.002 and row[0] <= 0.1) and (row[1] > 30000000 and row[1] <= 40000000)):
                row[2] = 11
                row[3] = 4.0
                row[4] = 0.05
                row[5] = 0.3
            elif ((row[0] > 0.002 and row[0] <= 0.1) and row[1] > 40000000):
                row[2] = 12
                row[3] = 4.5
                row[4] = 0.05
                row[5] = 0.35
            elif (row[0] > 0.01 and row[1] <= 1000000):
                row[2] = 13
                row[3] = 0.5
                row[4] = 0.1
                row[5] = 0.2
            elif (row[0] > 0.01 and (row[1] > 1000000 and row[1] <= 10000000)):
                row[2] = 14
                row[3] = 1.0
                row[4] = 0.1
                row[5] = 0.3
            elif (row[0] > 0.01 and (row[1] > 10000000 and row[1] <= 20000000)):
                row[2] = 15
                row[3] = 2.0
                row[4] = 0.1
                row[5] = 0.4
            elif (row[0] > 0.01 and (row[1] > 20000000 and row[1] <= 30000000)):
                row[2] = 16
                row[3] = 3.0
                row[4] = 0.1
                row[5] = 0.5
            elif (row[0] > 0.01 and (row[1] > 30000000 and row[1] <= 40000000)):
                row[2] = 17
                row[3] = 4.0
                row[4] = 0.1
                row[5] = 0.6
            elif (row[0] > 0.01 and row[1] > 40000000) :
                row[2] = 18
                row[3] = 4.5
                row[4] = 0.1
                row[5] = 0.7

            # Updating the cursor with the updated list
            cursor.updateRow(row)

# calculating soil depth based on min and max depth set at the begining          
# Org: Pacific Northwest National Laboratory
# ORIG-DATE:    Apr-2017
#Updated by: Madan Raj Bhandari (2024-5-09)
#This script create a average soil depth based on min, max depth and slope calculated from DEM. Area of high slope may have less soil depth and vice-versa
# mindepth - the minimum depth of the soil (this is a floor)
# maxdepth - the maximum depth of the soil (this will never be exceeded)
# wtslope - the relative weighting for the slope
# wtsource - the relative weighting for the source area
# wtelev - the relative weighting for the elevation
# maxslope - anything greater than this will create the slope function = 1
# maxsource - anything greater than this will create the source function = 1
# maxelev - anything greater than this will create the elevation function = 1
# powslope - raise the slope fraction by this power
# powsource - raise the source area fraction by this power
# powelev - raise the elevation fraction by this power

if arcpy.Exists(soildepth):
    print('Soil depth file already exists, deleting and creating a new one')
    arcpy.Delete_management(soildepth)
else:
    print('soil depth file not provided, creating a map')

#Starting calculation 
print ('starting soil depth calculation')
import arcpy
from arcpy import env
from arcpy.sa import *
import arcgisscripting
import os

# value can be modified or change based on individual scenario
wtslope = 0.6
print('The relative weighting for slope is '+ str(wtslope))
wtsource = 0.0
print('The relative weighting for source area is '+str(wtsource))
wtelev = 0.3
print('The relative weighting for elevation is '+ str(wtelev))
maxslope = 35.
maxsource = 10000000.
maxelev = 2500
powslope = .25
powsource = 1.
powelev = .75

print('All value read in')
totalwt = wtslope + wtsource + wtelev
print(str(totalwt))
     
if totalwt != 1.0:
    print('the weights must add up to 1.')
else:
    print('start creating file')
     
 
    slopegrid=Slope(dem,"DEGREE",1)
    slopegrid.save("slopegrid")
        
    tempsrc = Con(Raster(flowacc)>maxsource,maxsource,flowacc)
    tempelev = Con(Raster(dem)>maxelev,maxelev,dem)
    tempslope = Con(Raster("slopegrid")>maxslope,maxslope,slopegrid)

    print('calculate soil map')  
    soil_depth = mindepth + (maxdepth - mindepth)*\
                (wtslope * (1.0 - Power((tempslope/ maxslope), powslope))+ \
                 (wtsource * Power( (tempsrc/ maxsource), powsource))+\
                 wtelev *(1.0 - Power( (tempelev / maxelev), powelev)))
   
    soil_depth.save(soildepth)
        
    print('deleting temp files')  
    arcpy.Delete_management(tempsrc)
    arcpy.Delete_management(tempelev)       
    arcpy.Delete_management(tempslope)
    arcpy.Delete_management(slopegrid)


# Soil Detpth calculation for channel segment 

pol_ras="pol_ras"

if arcpy.Exists(pol_ras):
    print('Table already exists, deleting and creating a new one')
    arcpy.Delete_management(pol_ras)

arcpy.PolylineToRaster_conversion (streamnetwork, "arcid", pol_ras, "", "", env.cellSize)
soild_table="soildepth_zonal"

if arcpy.Exists(soild_table):
    arcpy.Delete_management(soild_table)
    print('zonal statistics table already exists, delete and create new')



ZonalStatisticsAsTable (pol_ras, "VALUE", soildepth, soild_table, "DATA", "MINIMUM") #calculates and generate soildepth for stream

arcpy.JoinField_management(streamnetwork,"arcid",soild_table ,"VALUE","MIN") #if this code doesnot work and if value are not presented in MIN column use manual techniques to insert value in MIN column

allnodes=path+"output.gdb/allnodes"    #creating allnodes feature
arcpy.FeatureVerticesToPoints_management(streamnetwork, allnodes, "ALL")
ExtractMultiValuesToPoints(allnodes, [[soildepth, "soildepth"]], "BILINEAR")

arcpy.Statistics_analysis(allnodes,"allnodes_statistics","soildepth MIN","arcid") #if error occur due to input fields containing null value use "summary statistics"geoprocessing tool to calculate statistical summary manually

arcpy.JoinField_management(streamnetwork,"arcid","allnodes_statistics","arcid","MIN_soildepth")

arcpy.AddField_management (streamnetwork, "segdepth", "FLOAT", 8, 2)
fields=['MIN_soildepth','MIN','segdepth']

with arcpy.da.UpdateCursor(streamnetwork, fields) as cursor:
    for row in cursor:
        if row[1] is None:
            row[2]=row[0]
        elif (row[0] < row[1]):
            row[2]=row[0]
        elif (row[0] > row[1]):
            row[2]=row[1]
        else:
            row[2]=row[0]
        cursor.updateRow(row)

arcpy.DeleteField_management(streamnetwork, "MIN")
arcpy.Delete_management(allnodes)
arcpy.Delete_management("allnodes_statistics")
#calculating effective depth
arcpy.CalculateField_management (streamnetwork, "effdepth", "0.95*float(!segdepth!)", "PYTHON3","")

arcpy.DeleteField_management(streamnetwork, "MIN_soildepth")

#Generating stream network file for running DHSVM this file will contain stream segment properties
    
if os.path.exists("stream.network.dat"):
    print('stream.network.dat already exists, delete and create new')
    os.remove("stream.network.dat")
    print('stream.network.dat sucessfully deleted')

print('creating new stream.network.dat file')


# Selecting fields for the stream network

arr_export=arcpy.da.TableToNumPyArray(streamnetwork,('arcid','segorder','slope','Shape_Length','chanclass','downarc'))

np.savetxt(path+'stream.network.dat', arr_export,fmt="%5d %3d %11.5f %17.5f %3d %7d")



#Creating stream map file for DHSVM
# This python script is created based on original AML scripts createstreamnetwork.aml as part of DHSVM. Reason for 4 neighbour: make the computation of surface and subsurface flow pathways consistent with the digital 
# elevation model networks (DEMON) described by Costa-Cabral and Burges (1994).

import arcpy
from arcpy import env
from arcpy.sa import *
import arcgisscripting
import os
import math
import numpy as np

rad2deg= 180.0 / math.pi

# Convert DEM to numpy array with no data value set to -9999
demarray = arcpy.RasterToNumPyArray(dem,nodata_to_value=-9999)
demraster=arcpy.sa.Raster(dem)

#Raster properties
elevcellx = arcpy.GetRasterProperties_management(dem, "CELLSIZEX")
deltax=float(elevcellx.getOutput(0))
elevcelly = arcpy.GetRasterProperties_management(dem, "CELLSIZEY")
deltay=float(elevcelly.getOutput(0))
    
lowerLeft = arcpy.Point(demraster.extent.XMin,demraster.extent.YMin)
    
# Convert flow accumulation to numpy array with no data value set to -9999
accarr = arcpy.RasterToNumPyArray(flowacc,nodata_to_value=-9999)
minloc=np.unravel_index(accarr.argmax(), accarr.shape)
accarr[minloc[0]][minloc[1]]=0
outlet=arcpy.NumPyArrayToRaster(accarr,lowerLeft,deltax,deltay)

size=demarray.shape  # size=[Number of Rows, Number of Columns]

# Extend boundary of dem input by add empty row/column before/after first/last row/column
demarr=np.zeros((size[0]+2,size[1]+2))
slopearr=np.zeros((size[0],size[1]))
flag=np.zeros((size[0],size[1]))
tmp=np.zeros((size[0],size[1]))

# Change array value to nodata value in DEM array
demarr[:][:]=-9999
for ii in range(1,size[0]+1):
    demarr[ii][1:(size[1]+1)]=demarray[ii-1][:]

# Calculate slope for extended DEM array
for i in range(1,size[0]+1):  # loop through i-row
    for j in range(1,size[1]+1): # loop through j-col
        if (demarr[i][j+1]==-9999) and (demarr[i][j-1]==-9999):            #calculating dx
            dzdx=0.0
        elif demarr[i][j-1]==-9999: #aspect(0-180) or(0,pi),dzdz>=0
            dzdx=max((demarr[i][j]-demarr[i][j+1])/(deltax),0.0)
        elif demarr[i][j+1]==-9999: #aspect(180-360) or(-pi,0),dzdz<=0
            dzdx=min((demarr[i][j-1]-demarr[i][j])/(deltax),0.0)
        else:
            dzdx=(demarr[i][j+1]-demarr[i][j-1])/(2*deltax)

# calculating dy
        if (demarr[i+1][j]==-9999) and (demarr[i-1][j]==-9999):
            dzdy=0.0
        elif demarr[i-1][j]==-9999: #aspec (pi/2,pi)||((-pi,-pi/2)), dzdy<=0
            dzdy=min((demarr[i+1][j]-demarr[i][j])/(deltay),0.0)
        elif demarr[i+1][j]==-9999: #aspec (0, pi/2)||((-pi/2,0)), dzdy>=0
            dzdy=max((demarr[i][j]-demarr[i-1][j])/(deltay),0.0)
        else:
            dzdy=(demarr[i+1][j]-demarr[i-1][j])/(2*deltay)
                
        slopearr[i-1][j-1]=math.sqrt(math.pow(dzdx,2)+math.pow(dzdy,2))
        if ((dzdx==0) and (dzdy==0)):
            flag[i-1][j-1]=1        
        else:
            flag[i-1][j-1]=0
            tmp[i-1][j-1]=int(math.atan2(dzdx,dzdy) * rad2deg +360)

flag
tmp
FlagRaster = arcpy.NumPyArrayToRaster(flag,lowerLeft,deltax,deltay)
TmpRaster = arcpy.NumPyArrayToRaster(tmp,lowerLeft,deltax,deltay)

tmp2=Aspect(dem) 
junk=Con(FlagRaster==1,tmp2,TmpRaster)
tmpaspect=Con(junk>360,junk%360,junk)

tmpslope = arcpy.NumPyArrayToRaster(slopearr,lowerLeft,deltax,deltay)
sloperas=SetNull(dem==-9999,tmpslope)
aspectras=SetNull(dem==-9999,tmpaspect)
slope="gridslope_4d"
aspect="gridaspect" 
sloperas.save(slope)
aspectras.save(aspect) 
    


#Creating row column polygon feature for stream map 
import arcpy
from arcpy import env
from arcpy.sa import *
import arcgisscripting
import math
import numpy as np

env.mask = MaskExt

env.cellSize=dem
arcpy.env.extent=dem

demraster=arcpy.sa.Raster(dem)
RightRaster=CreateConstantRaster (1, "INTEGER", env.cellSize, demraster.extent) 
UpRaster=CreateConstantRaster (4, "INTEGER", env.cellSize, demraster.extent)
    
col=FlowAccumulation(RightRaster)
col.save("colraster")
    
rownum=FlowAccumulation(UpRaster)
rownum.save("rowraster")

print('row and column completed')

# Create Grid-Points and Grid-Polygon
if arcpy.Exists("cb"):
    arcpy.Delete_management("cb")
        
print('joining row map and col map')
rowcolcomb=Combine(["colraster","rowraster"])
rowcolcomb.save("cb")

print('rowcolcomb done')

if arcpy.Exists("rowcolpoly.shp"):
    arcpy.Delete_management("rowcolpoly.shp")
#converting Raster to polygon
print('converting raster to polygon, may take couple of minutes; sometime system may crash due to large file size')
arcpy.RasterToPolygon_conversion(rowcolcomb, "rowcolpoly.shp", "NO_SIMPLIFY", "VALUE")

print('Joinning table, this may take a while...')

joinfields = ['VALUE', 'COLRASTER', 'ROWRASTER']
joindict = {}
with arcpy.da.SearchCursor("cb", joinfields) as rows:
    for row in rows:
        joinval = row[0]
        val1 = row[1]
        val2 = row[2]
        joindict[joinval]=[val1, val2]
del row, rows

print('adding fields to polygon feature (i.e. rowcolpoly.shp)')
arcpy.AddField_management("rowcolpoly.shp", "ROW", "LONG")
arcpy.AddField_management("rowcolpoly.shp", "COL", "LONG")

  
targfields = ['GRIDCODE', 'COL', 'ROW']
with arcpy.da.UpdateCursor("rowcolpoly.shp", targfields) as recs:
    for rec in recs:
        keyval = rec[0]

        if keyval in joindict:
            rec[1] = joindict[keyval][0]
            rec[2] = joindict[keyval][1]
        else:
            rec[1] = 0
            rec[2] = 0
        recs.updateRow(rec)
del rec, recs
    
print('Join table finished')
    
print('deleting unwanted file')    
arcpy.Delete_management(col)
arcpy.Delete_management(rownum)
arcpy.Delete_management(rowcolcomb)
arcpy.Delete_management(RightRaster)
arcpy.Delete_management(UpRaster)

#creating a outcover file for stream map and computing geometric interection 
if arcpy.Exists('/output.gdb/outcover'):
    arcpy.Delete_management('/output.gdb/outcover')

outcover=".//output.gdb//outcover"
arcpy.Intersect_analysis([streamnetwork,"rowcolpoly.shp"], outcover, "", "", "")


#roadaspects calculation
#Part of python version createstreamnetwork. Computes the aspect of stream within each grid cell 
import arcpy
import math

arcpy.AddField_management (outcover, "rd_aspect", "FLOAT", 12, 5)      #segment Aspect for stream map
arcpy.AddField_management (outcover, "rd_st_len", "FLOAT", 12, 5)      #used for road Network (if applicable)
arcpy.AddField_management (outcover, "rd_efflen", "FLOAT", 12, 5)      #used for road network (if applicable)

#calculating aspect (degree) of stream segment
print('calculating stream segment aspect')

expression1 = "GetGeographicalDegrees( !SHAPE! )"
    
codeblock1 = """
def GetGeographicalDegrees(shape):
    import math
    radian = math.atan2(shape.firstpoint.y - shape.lastpoint.y, shape.firstpoint.x - shape.lastpoint.x)
    rad2deg = 180.0 / math.pi
    degrees = round((radian * rad2deg + 360) % 360)
    return degrees
"""

arcpy.CalculateField_management (outcover, "rd_aspect", expression1, "PYTHON3", codeblock1)


expression2 = "CalcStLen( !SHAPE!)"
codeblock2 = """
def CalcStLen(shape):
    import math
    length = math.sqrt((shape.lastpoint.y - shape.firstpoint.y) ** 2 +
                       (shape.lastpoint.x - shape.firstpoint.x) ** 2)
    return length
"""

arcpy.CalculateField_management (outcover, "rd_st_len", expression2, "PYTHON3", codeblock2)

expression3 = "CalcEffectLen( !SHAPE!,!rd_st_len! ,!rd_aspect! )"
codeblock3 = """
def CalcEffectLen(shape, length, rd_aspect):
    import math
    eff_len = math.fabs(length * math.sin(rd_aspect * math.pi / 180))
    return eff_len
"""

arcpy.CalculateField_management (outcover, "rd_efflen", expression3, "PYTHON3", codeblock3)

#Generating stream Map file (output)
import arcpy
import numpy as np
import os
import datetime

print('generating stream map file')
if os.path.exists("stream.map.dat"):
    os.remove("stream.map.dat")
    print('stream.map.dat sucessfully deleted')

print('creating new stream.map.dat file')

info = '###### This file has been automatically generated #####'
info += '\n######             EDIT WITH CARE!!!              #####'
info += '\n# Generated: '+str(datetime.datetime.now())
info += '\n# Created by streammapfile.py'
info += '\n# Workspace '+str(path)
info += '\n#                   Segment  Cut/Bank     Cut     Segment'
info += '\n#  Col  Row  ID      Length   Height     Width     Aspect   SINK?'
info += '\n#                     (m)      (m)        (m)       (d)    (optional)'
info += '\n# \n'

arr_export=arcpy.da.TableToNumPyArray(outcover,('COL','ROW','arcid','Shape_Length','effdepth','effwidth','rd_aspect'))

with open(path+'stream.map.dat', 'w') as f:
    f.write(info)
    np.savetxt(f, arr_export,fmt="%5d %5d %5d %11.4f %10.4f %9.4f %10.4f")

#cleaning the unwanted files
print('deleting unwanted files and completing the run')
arcpy.Delete_management(flowdir)
arcpy.Delete_management(flowacc)
arcpy.Delete_management(slope)
arcpy.Delete_management(aspect)
arcpy.Delete_management(tmpelev)
arcpy.Delete_management(MaskExt)
arcpy.Delete_management(local)
arcpy.Delete_management(soild_table)
arcpy.Delete_management(outcover)
rcpy.Delete_management(nodestart)
arcpy.Delete_management(nodeend)
arcpy.Delete_management(tmpacc)
arcpy.Delete_management(Pol_ras)
