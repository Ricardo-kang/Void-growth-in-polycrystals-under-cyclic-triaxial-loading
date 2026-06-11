#导入python模块
import os
from ovito.io import import_file
from ovito.data import SurfaceMesh
from ovito.data import SimulationCell
from ovito.modifiers import ConstructSurfaceModifier
import math
import time


#更改路径
path_input=os.getcwd()
#path_input=r'E:\0 MD-code\0 Polycrystals-Cu\240-27-20-ori-3\20221019-void-240-grain=27-R=20-T=3-model-6-ori-3\result'
filename_input='tension.*.dat'

path_output=os.getcwd()
#path_output=r'E:\0 MD-code\0 Polycrystals-Cu\240-27-20-ori-3\20221019-void-240-grain=27-R=20-T=3-model-6-ori-3\1 date'
filename1_output='Volume-fraction.txt'

#读入lammps的轨迹文件，导入pipeline(ovito计算流), 
pipeline = import_file(path_input+'/'+filename_input)

f = open(path_output+'/'+filename1_output,'w')
#f.write("Step Strain Volume_fraction_Void\n")
f.close()

#逐个读取轨迹文件，并进行DXA分析与输出
for frame in range(pipeline.source.num_frames):
        time_start = time.perf_counter()
        print('----------------------------------------------------------------------------------------------') 
        data = pipeline.compute(frame)
        pipeline.modifiers.clear()
        #当前文件时间步:
        step=data.attributes['Timestep']
        print('(1)当前时间步： '+str(step)+'\n')
        
        #计算当前应变x：
        lx = data.cell[0,0]
        print('(2)当前模拟盒x向尺寸：'+str(lx)+' Am\n')
        print(str(lx))
        if frame == 0:
           lx0 = lx
        strainx = (lx-lx0)/lx0
        print('(3)当前应变：         '+str(strainx)+'\n')


        #计算当前应变y：
        ly = data.cell[1,1]
        print('(2)当前模拟盒x向尺寸：'+str(ly)+' Am\n')
        print(str(ly))
        if frame == 0:
           ly0 = ly
        strainy = (ly-ly0)/ly0
        print('(3)当前应变：         '+str(strainy)+'\n')


        #计算当前应变z：
        lz = data.cell[2,2]
        print('(2)当前模拟盒x向尺寸：'+str(lz)+' Am\n')
        print(str(lz))
        if frame == 0:
           lz0 = lz
        strainz = (lz-lz0)/lz0
        print('(3)当前应变：         '+str(strainz)+'\n')

        straineq = math.sqrt(((strainx - strainy)**2+(strainy - strainz)**2+(strainz - strainx)**2)*2/9)


        #定义PTM模块,设置radius值
        pipeline.modifiers.append(ConstructSurfaceModifier(
            method = ConstructSurfaceModifier.Method.AlphaShape,
            radius = 2.9,
            identify_regions = True))
        #pipeline.modifiers.append(ConstructSurfaceModifier(radius = 2.9))
        #pipeline计算启动
        data = pipeline.compute(frame)
        #print("Solid volume: %f" % data.attributes['ConstructSurfaceMesh.filled_volume'])
        vol = data.attributes['ConstructSurfaceMesh.cell_volume']
        Solid_volume = data.attributes['ConstructSurfaceMesh.filled_volume']
        print( '(4)Cell volume： '+str(vol)+' Am^3\n')
        print('(6)Solid volume:  ' +str(Solid_volume)+' Am^3\n')

        Void_volume = vol-Solid_volume
        print('(7)Void volume:  ' +str(Void_volume)+' Am^3\n')

        Solid_volume_fraction = Solid_volume / vol
        Void_volume_fraction = 1-Solid_volume_fraction
        print("(8)Void volume fraction: %f" % Void_volume_fraction)

        #输出并关闭已打开的修正:
        f = open(path_output+'/'+filename1_output,'a')
        f.write(str(step)+' '+str(strainx)+' '+str(strainy)+' '+str(strainz)+' '+str(straineq)+' '+str(vol)+' '+str(Solid_volume)+' '+str(Void_volume_fraction)+' '+'\n')
        f.close()
        del pipeline.modifiers[0]

        time_end = time.perf_counter()
        print("用时 : {} s".format(time_end - time_start))
