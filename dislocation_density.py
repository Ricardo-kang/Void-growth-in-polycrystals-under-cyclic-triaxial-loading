#导入python模块
import os
from ovito.io import import_file
from ovito.modifiers import DislocationAnalysisModifier
from ovito.modifiers import SelectTypeModifier
from ovito.modifiers import ConstructSurfaceModifier
from ovito.data import DislocationNetwork
from ovito.data import SimulationCell
import math
import time

#更改路径
path_input=os.getcwd()
#path_input=r'E:\0 MD-code\0 Polycrystals-Cu\240-27-20-ori-3\20221019-void-240-grain=27-R=20-T=3-model-6-ori-3\result'
filename_input='tension.*.dat'

path_output=os.getcwd()
#path_output=r'E:\0 MD-code\0 Polycrystals-Cu\240-27-20-ori-3\20221019-void-240-grain=27-R=20-T=3-model-6-ori-3\1 date'
filename1_output='dislocation_density_sc_2.txt'

#读入lammps的轨迹文件，导入pipeline(ovito计算流), 
pipeline = import_file(path_input+'/'+filename_input)

f = open(path_output+'/'+filename1_output,'w')
#f.write("Strain-TDD-SDD-PDD-SrDD-HDD <Units:A^-2>\n")
f.close()

#逐个读取轨迹文件，并进行DXA分析与输出
for frame in range(pipeline.source.num_frames):
        time_start = time.perf_counter()
        print('----------------------------------------------------------------------------------------------') 
        data = pipeline.compute(frame)

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

        #（1）分析双层材料中总的位错密度
        #计算位错长度：
        modifier1 = DislocationAnalysisModifier()
        modifier1.input_crystal_structure = DislocationAnalysisModifier.Lattice.FCC       # 晶格类型
        pipeline.modifiers.append(modifier1)
        data = pipeline.compute(frame)
        Tll  = data.attributes['DislocationAnalysis.total_line_length']   # total
        #Oll  = data.attributes['DislocationAnalysis.length.other']      #other
        Sll  = data.attributes['DislocationAnalysis.length.1/6<112>']     # Shockley_dislocation
        Pll  = data.attributes['DislocationAnalysis.length.1/2<110>']     # Perfect_dislocation
        Srll = data.attributes['DislocationAnalysis.length.1/6<110>']     # Stair_rod_dislocation
        Hll  = data.attributes['DislocationAnalysis.length.1/3<100>']     # Hirth_dislocation
        #Fll  = data.attributes['DislocationAnalysis.length.1/3<111>']     #Frank

        print('(A1)双层材料各类位错长度：(Am)\n')
        print('    总位错长度：      '+str(Tll)+'\n')
        #print('    other位错长度：      '+str(Oll)+'\n')
        print('    Shockley位错长度：'+str(Sll)+'\n')
        print('    理想位错长度：    '+str(Pll)+'\n')
        print('    压杆位错长度：    '+str(Srll)+'\n')
        print('    Hirth位错长度：   '+str(Hll)+'\n')
        #print('    Frank位错长度：   '+str(Fll)+'\n')
       
        #计算体积：
        
        vol = data.cell.volume
        print( '(A2)双层材料的体积： '+str(vol)+' Am^3\n')

        #计算位错密度：
        TDD  = Tll * 1000 / vol
        #ODD  = Oll * 1000 / vol
        SDD  = Sll * 1000 / vol
        PDD  = Pll * 1000 / vol
        SrDD = Srll * 1000/ vol
        HDD  = Hll * 1000/ vol
        #FDD  = Fll * 1000/ vol
        print('(A3)双层材料中各类位错密度：(Am^-2)\n')
        print('    总位错密度：      '+str(TDD)+'\n')
        #print('    Other位错密度：      '+str(ODD)+'\n')
        print('    Shockley位错密度：'+str(SDD)+'\n')
        print('    理想位错密度：    '+str(PDD)+'\n')
        print('    压杆位错密度：    '+str(SrDD)+'\n')
        print('    Hirth位错密度：   '+str(HDD)+'\n')
        #print('    Frank位错密度：      '+str(FDD)+'\n')

        
        #输出并关闭已打开的修正:
        f = open(path_output+'/'+filename1_output,'a')
        f.write(str(step)+' '+str(strainx)+' '+str(strainy)+' '+str(strainz)+' '+str(straineq)+' '+str(TDD)+' '+str(SDD)+' '+str(PDD)+' '+str(SrDD)+' '+str(HDD)+' '+'\n')
        f.close()
        del pipeline.modifiers[0]

        time_end = time.perf_counter()
        print("用时 : {} s".format(time_end - time_start))
