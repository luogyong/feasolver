LGYWORKS

iFsolver输入文件

帮助文档(V0.1-2022-08-09)

# 介绍

## 主要功能

iFsolver的主要功能是进行有限元问题的求解，其输入文件为sinp文件。求解完成后输出后处理tecplot文件。

## 代码

<https://github.com/luogyong/feasolver>

## 输入的基本格式

关键词+参数

并非每个关键词都要输入，有需要输入。除非特别说明，各关键词的顺序也无关紧要。

# iFsolver输入文件各关键词输入格式说明

## Title

表 1 title输入格式说明表

|  |  |  |
| --- | --- | --- |
| **关键词** | **TITLE** | |
| 功能 | 模型的名称 | |
| 输入格式： | | |
| Title  *string* | | |
| 示例： | | |
| Title  aniso\_1D\_flow | | |
| 字段说明： | | |
| *string* | | 字符串，模型的名称 |

## Node

表 2 Node输入格式说明表

|  |  |  |
| --- | --- | --- |
| **关键词** | **Node** | |
| 功能 | 输入模型的节点坐标 | |
| 输入格式： | | |
| Node,*NUM*=I,*DATAPACKING*=0|1,*DIMENSION*=2|3,*isporeflow*=0|1  *X*,*Y*[,*Z*] | | |
| 示例： | | |
| NODE,NUM=16, DATAPACKING=1, DIMENSION=2  0.1465762786596120 0.2506144708076132  0.2500000000000000 0.0000000000000000  …… | | |
| 字段说明： | | |
| *NUM* | | 节点数 |
| *DATAPACKING* | | 输入格式  =1(默认)，按点的顺序输入{x1,y1,z1},{x2,y2,z2},…  =0，按坐标的顺序输入{x1,x2,…},{y1,y2,…},{z1,z3,…} |
| *isporeflow* | | 是否为孔隙管流计算  =0，否（默认）  =1，是 |
| *x,y,z* | | 点的坐标  当dimension=2时，二维模型，输入x,y  当dimension=3时，三维模型，输入x,y,z |

## ELEMENT

表 3 ELEMENT输入格式说明表

|  |  |  |
| --- | --- | --- |
| **关键词** | **ELEMENT** | |
| 功能 | 输入单元信息 | |
| 输入格式： | | |
| ELEMENT,*num*=I,*set*=I,*et*=C,*matid*=I,coupleset=I3,*sf*=I,*istopo*=I,*title*=C  *n*1,*n*2,…*nn*,[*property*] //各单元的节点编号 | | |
| 示例： | | |
| ELEMENT,*num*=18, *set*=201, *et*=cpe3\_spg, *matid*=1, *coupleset*=201, *sf*=0, *istopo*= 0, *title*=model  1 2 3  4 1 5  3 2 6  …… | | |
| 字段说明： | | |
| *num* | | 单元数 |
| *set* | | 单元组编号 |
| *et* | | 单元类型  目前支持以下单元类型：  渗流单元：  二维：CPE3\_SPG, CPE6\_SPG,CPE15 \_SPG  轴对称：CAX3\_SPG, CAX6\_SPG, CAX15\_SPG  三维：TET4\_SPG,TET10\_SPG, PRM6\_SPG, PRM15\_SPG  无厚度：ZT4\_SPG2, ZT4\_SPG, ZT6\_SPG2, ZT6\_SPG  井流单元：WELLBORE, WELLBORE\_SPG,  球状流单元:sphflow, semi\_sphflow  管流单元：PIPE2,poreflow  一般力学：  平面应变单元：CPE3,CPE6, CPE15,  平面应力单元：CPS3,CPS6,CPS15  三维单元：PRM6,PRM15,TET4,TET10  轴对称单元：CAX3,CAX6,,CAX15  结构单元：  杆单元：BAR，BAR2D  梁单元：BEAM，BEAM2D  壳单元：SHELL3  弹簧单元：SPRINGX, SPRINGY, SPRINGZ, SPRINGMX,  SPRINGMY, SPRINGMZ  土土弹簧单元：soilspringx, soilspringy, soilspringz |
| *matid* | | 单元材料号 |
| *coupleset* | | 与此单元组耦合的单元组号，默认自身，即没有耦合单元组 |
| *sf* | | 步函数号，默认0。 |
| *istopo* | | 大多数情况下为可输入0。  当单元为wellbore单元，ISTOPO=1时表示由高次线单元分解形成的一次的wellbore单元，（当井周单元为线性单元时，不存在此问题，*istopo*=0）。因后面单元拓扑邻接分析时，只分析高次单元的边(即忽略中间节点)，导致由含内部节点的单元的边不体现，为解决此问题，输出对应高次单元的端节点，利用此端节点进行拓扑分析。略显麻烦。  当*istopo=0,*要通过*property*输入*wellbore*对应高次单元的端节点。 |
| *title* | | 单元组名称 |
| *n1,n2,..,nn* | | 单元节点号 |
| *property* | | 可选参数，存在下列情况下输入：   1. 当单元为shell3时，输入单元的厚度 2. 当单元为semi\_sphflow时，输入半球流的方向矢量 3. 当单元为zt6\_spg, zt4\_spg, zt6\_spg2, zt4\_spg2时，可输入各单元的材料号 4. 当为poreflow/pipe2时，可以输入喉径(**不是水力半径**)、单元阻力(单元渗透系数的倒数) |
| 注意 | | 此关键词可重复出现 |

## MATERIAL

表 4 MATERIAL输入格式说明表

|  |  |  |
| --- | --- | --- |
| **关键词** | **MATERIAL** | |
| 功能 | 输入材料参数 | |
| 输入格式： | | |
| MATERIAL,*matid*=I, *name*=Char[*type*=Char]  *property* | | |
| 示例： | | |
| MATERIAL,*matid*=1,name=soil  1,1.0,1 | | |
| 字段说明： | | |
| *matid* | | 材料号 |
| *name* | | 材料名 |
| *type* | | 材料类型，可选参数。可输入的字符为：  Elastic：当单元为固体/结构单元时的默认值。目前可处理的非线性材料模型有：   1. Mises：Mises弹理想塑性材料   2) MC：摩尔库伦材料  3) eip\_bar：弹理想塑杆  4) eip\_beam：弹理想塑梁  linear\_spg：当单元为渗流单元时的默认值(饱和土渗流模型)。渗流模型还有：1)step\_spg ,饱和土渗流模型;  2)vg\_spg, van genuchten非饱和土渗流模型；  3) lr\_spg, Leong and Rahardjo非饱和土渗流模型；  4) exp\_spg, 指数非饱和土渗流模型 |
| *property* | | 材料参数。根据材料和单元不同按下列顺序输入：  1)当为弹性材料时：  弹性模量E, 泊松比v  2)当为渗流材料时：  k1,k2[,k3][,AnisoTransM,0,0, alpha,n,m,Mv,sita\_s,sita\_r,rw,thickness]  其中，  k1,k2,k3为三个主渗透方向的参数；  AnisoTransM，主渗流方向与整体坐标方向不一致时，其转化矩阵(通过关键词coordinate输入)的编号；  alpha,n,m为van genuchten和Leong and Rahardjo非饱和土渗流模型；  Mv，压缩系数  sita\_s， 饱和体积含水量  sita\_r，残余体积含水量  rw，水的重度  thickness ,无厚度单元的厚度默认为1.  注意：zt4\_spg and zt6\_spg, 实际为一维流单元，只用k(ndimension)表示其垂直于单元边或面的渗透系数，thickness为防渗墙厚度。  3) Mises材料  弹性模量E, 泊松比v,单轴抗拉强度sigma\_y，土的重读rw;  注意：在平面应变模型中， sigma\_y=3\*Cu (Cu, undrained shear strength ) ；在轴对称模型中，sigma\_y=2\*Cu.   1. MC材料   弹性模量E, 泊松比v, 黏聚力C, 摩擦角phi; 剪胀角D；土重rw;   1. 当为Bar单元材料时   E弹模,截面积A, 梁截面y'(局部坐标)的高度hy, 梁截面z'(局部坐标)的高度hz, 所能承受的最小的轴向压力MinN, 所能承受的最大的轴向压力MinN MaxN   1. 当为beam单元材料时   E弹模,截面积A, 泊松比u,J,Iz,Iy, hy, hz, MinN, MaxN,MinMx,MaxMx, MinMy,MaxMy, MinMz,MaxMz(分别梁局部坐标下的为轴力和弯矩的极限值)   1. 当材料为管(井)流单元时 2. R(井半径,注意不是水力半径) 3. Temp(温度), 4. Kr (relative roughness of the inner surface of the pipe ) 5. PIPE-FLOW MODEL:   =0,Darcy(no porous effect,DEFAULT);  =1,Siwon;  =2,OUYang EFFECT;  =3,Input by user;  =4,laminar flow, no need to update the friction factor.   1. 泥皮的厚度和渗透系数的比L/K.=0(默认,不考虑井损)). 2. f(Darcy Friction factor(井壁摩阻系数),=0,由计算定(默认)，一般为0.02-0.03左右,当model=3时输入。). 3. POROSITY OF THE WELLBORE(当model=1时输入). |
| 注意 | | 此关键词可重复出现 |

## BC

表 5 BC输入格式说明表

|  |  |  |
| --- | --- | --- |
| **关键词** | **BC** | |
| 功能 | 输入模型边界条件 | |
| 输入格式： | | |
| BC,*num*=I[,*isinc*=I,*stepfunc*=I,*iswellhead*=I,*spg\_isdual*=I, *ssp\_onepile*=I] //命令行  *node,dof,value [,stepfunc,spg\_isdual,ssp\_onepile,isinc]* //数据行 | | |
| 示例： | | |
| BC,NUM=6,ISINC=0  16 4 1.5000000000000000 0 0  14 4 1.5000000000000000 0 0  12 4 1.5000000000000000 0 0  5 4 0.5000000000000000 0 0  4 4 0.5000000000000000 0 0 | | |
| 字段说明： | | |
| *num* | | 此边界数据行数 |
| *isinc* | | 边界值是否为增量(多步计算时适用)  =0,value值当前步的边界值(默认)，  =1，value值为增量，当前步的边界值为前面各步值的和。 |
| *stepfunc* | | 步函数号，默认为0 |
| *iswellhead* | | 是否为自流井边界节点  0（默认），不是  1，是，表示此节点为减压自流井边界节点，节点流量只出不进，如果流入流量，则令此边界失效 |
| *spg\_isdual* | | 如果*isdual*=i(>0),则表示此自由度可能与出溢边界*Nseep*(i)重复，如果边界水头小于位置水头，则变为出溢边界。  默认为0. |
| *ssp\_onepile* | | 只对SSP单元节点起作用，标示这个作用是否是作用在其中一根钢板桩上,而不两根都作用。  /=0,作用在单根钢板桩上；  =0都作用（默认） |
| *node* | | 边界作用的节点号 |
| *dof* | | 节点自由度编号 |
| *value* | | 边界值  如果stepfunc>0,isinc=0时第istep步的边界值  *Vi*= value \* sf(stepfunc).factor(istep)  如果stepfunc>0,isinc=0时第istep步的边界值  *Vi*= Sum(value \* sf(stepfunc).factor(istep)) |
| 注意 | | 1. 当命令行与数据行出现相同的参数时，以数据行为准。 2. 此命令行可以重复出现 3. 对于水头边界(dof=4)，如果边界值小于其相应的位置水头，则认为该边界无效，不起作用。 4. 多步计算时，输入时要求每一步的水头均是总量，而不是增量. |

## LOAD

此关键词与BC的格式类似。

表 6 LOAD输入格式说明表

|  |  |  |
| --- | --- | --- |
| **关键词** | **LOAD** | |
| 功能 | 输入模型荷载信息 | |
| 输入格式： | | |
| Load,*num*=I[, *isinc*=I, *stepfunc*=I, *ssp\_onepile*=I] //命令行  *node, dof, value [,stepfunc, ssp\_onepile, isinc]* //数据行 | | |
| 示例： | | |
| Load,NUM=6,ISINC=0  16 4 1.5000000000000000 0 0  14 4 1.5000000000000000 0 0  12 4 1.5000000000000000 0 0  5 4 0.5000000000000000 0 0  4 4 0.5000000000000000 0 0 | | |
| 字段说明： | | |
| *num* | | 此荷载数据行数 |
| *isinc* | | 边界值是否为增量(多步计算时适用)  =0,value值当前步的边界值(默认)，  =1，value值为增量，当前步的边界值为前面各步值的和。 |
| *stepfunc* | | 步函数号，默认为0 |
| *ssp\_onepile* | | 只对SSP单元节点起作用，标示这个作用是否是作用在其中一根钢板桩上,而不两根都作用。  /=0,作用在单根钢板桩上；  =0都作用（默认） |
| *node* | | 边界作用的节点号 |
| *dof* | | 作用的自由度编号   1. XDIS 2. YDIS 3. ZDIS 4. H 5. SX(X转角) 6. SY 7. SZ |
| *value* | | 边界值  如果stepfunc>0,isinc=0时第istep步的边界值  *Vi*= value \* sf(stepfunc).factor(istep)  如果stepfunc>0,isinc=1时第istep步的边界值  *Vi*= Sum(value \* sf(stepfunc).factor(istep)) |
| 注意 | | 1. 当命令行与数据行出现相同的参数时，以数据行为准。 2. 此命令行可以重复出现 |

## SEEPAGE FACE

表 7 SEEPAGE FACE输入格式说明表

|  |  |  |
| --- | --- | --- |
| **关键词** | **SEEPAGE FACE** | |
| 功能 | 出溢面节点 | |
| 输入格式： | | |
| seepage face, *num*=I [,*sf*=I,*iswellbore*=I]  *n*1,*n*2,…*nnum* | | |
| 示例： | | |
| seepage, num=2, sf=0, iswellbore=0  1,2 | | |
| 字段说明： | | |
| *num* | | 节点数 |
| *sf* | | 步函数 |
| *iswellbore* | | iswellhead>0, 表示此节点为井壁出溢面，且此出溢点的流量，计入编号为iswellhead的井点流量。  =0，默认。 |
| *n*1,*n*2,…*nnum* | | 出溢面上的节点号 |
| 注意 | | 此关键词可重复出现 |

## Initial Value

此关键词与BC的格式类似。

表 6 Initial Value输入格式说明表

|  |  |  |
| --- | --- | --- |
| **关键词** | **IV** | |
| 功能 | 输入模型初始条件信息 | |
| 输入格式： | | |
| IV,*num*=I//命令行  *node, dof, value* //数据行 | | |
| 示例： | | |
| Load,NUM=6  16 4 1.5000000000000000  14 4 1.5000000000000000  12 4 1.5000000000000000  5 4 0.5000000000000000  4 4 0.5000000000000000 | | |
| 字段说明： | | |
| *num* | | 此初始条件数据行数 |
| *stepfunc* | | 步函数号，默认为0 |
| *node* | | 作用的节点号 |
| *dof* | | 作用的自由度编号   1. XDIS 2. YDIS 3. ZDIS 4. H 5. SX(X转角) 6. SY 7. SZ |
| *value* | | 初始值 |
| 注意 | | 1. 没有输入的节点，默认初始值为0。 2. 当进行孔隙网络淤堵模型分析时，可以将此命令改为“IC”,以输入各节点的初始体积浓度。此时，dof可以输入任意值。 |

## Solver

表 8 Solver输入格式说明表

|  |  |  |
| --- | --- | --- |
| **关键词** | | **SOLVER** |
| 功能 | | 求解控制参数 |
| 输入格式： | | |
| Solver, | | |
| 示例： | | |
| solver,*solver*=N\_R,*bfgm*=continuum,*niteration*=1000,*isfc*=1,*ftol*=1e-3， *wellmethod*=2,*wellaniso*=2 | | |
| 字段说明： | | |
| *solver* | 方程求解方法，默认为Newton-Raphson方法 | |
| *bfgm* | 应力积分方法，目前默认参数为continuum | |
| *niteration* | 最大的迭代计算次数 | |
| *isfc* | 以力/流量的平衡判断收敛性(默认) | |
| *ftol* | 力/流量的收敛准则 | |
| *i2n* | 积分点应力向节点应力的映射方法  =expo，默认；  =spr, 由节点周边的积分点进行插值  =avg, 由节点周边的积分点进行平均  =sht, 由最近的积分点平移 | |
| *i2nw* | 权，当i2n=avg时，  = WEIGHT\_ANGLE,以单元的角度为权  =weight\_nelt，默认，按单元的个数平均 | |
| *nopopup* | 是否采用静默的运行方式程序，即程序无错运行结束后退出。  =No，默认； | |
| *wellmethod* | 计算解析井流量采样方法。  =0，取样点为单元节点；  =1，为单元形心;  =2,在井线单元周边均布3层采样点(两端及中间),每排采样点数为nspwell，而球状流仍为单元节点。  =3,按球均匀采样（默认）  =4，在井线单元周边均布3层采样点(两端及中间),每排采样点数为nspwell，而球状流的方法为解析解（虚拟旁路单元） | |
| *wellaniso* | 水平各向异性井流解析流量的计算方法。  =0，directional K method；(默认)  =1, Charles R. Fitts method.(转化为各向同性材料进行)；  =2，王建荣方法，根据其公式进行计算 | |
| *isParasys* | 是否为参数敏感分析  =0（默认），不是  >0,是 | |
| *disf\_scale* | 位移放大倍数（tecplot显示变形后的图形用）  默认为1 | |
| *slope\_kscale* | 流线法边坡稳定分析参数  =2.0，默认。 | |
| *slope\_kbase* | 流线法参数  如果slope\_kbase<0,则kb取相对值，kb=abs(slope\_kbase)\*maxsfr;  否则，kb取绝对值，kb=slope\_kbase  =-0.2，默认。 | |
| *slope\_kratio* | 流线法参数  if(slope\_kratio>0),ky=kx/slope\_kratio,else,ky=input value. | |
| *slope\_mko* | 流线法参数  slope\_mko>0,表sfr要减掉初始ko应力场的sfrko.  假定ko=v/1-v,sxx=ko\*syy,szz=sxx,txy=txz=tyz=0 | |
| *slidedirection* | 边坡的滑移方向  =right,默认  =left | |
| *time\_unit* | 模型时间单位   1. 0, day 默认 2. 1, sec 3. 2, min 4. 3, hour | |
| *len\_unit* | 模型长度单位   1. 0, m默认 2. 1, dm 3. 2, cm 4. 3,mm 5. 4,km | |
| *well\_bottom\_type* | 井底的形状   1. 0, 平底（默认） 2. 非0，球面状 | |

## OUTVAR

表 9 OUTVAR输入格式说明表

|  |  |  |
| --- | --- | --- |
| **关键词** | **OUTVAR** | |
| 功能 | 控制计算输出量 | |
| 输入格式： | | |
| OUTVAR,[*spg,dis,disf,stress,mises,strain,pstrain,eeq,peeq,moment,rotate,force,sfr,vsf*] | | |
| 示例： | | |
| outvar,spg | | |
| 字段说明： | | |
| *spg* | | 输出渗流计算相关的量 |
| *dis* | | 输出节点位移 |
| *disf* | | 输出位移后节点坐标 |
| *stress* | | 输出节点应力 |
| *strain* | | 输出节点应变 |
| *pstrain* | | 输出节点塑性应变 |
| *mises* | | 输出节点miese应力 |
| *eeq* | | 输出节点等效应变 |
| *peeq* | | 输出节点等效塑性应变 |
| *moment* | | 输出节点弯矩 |
| *rotate* | | 输出节点转角 |
| *force* | | 输出节点力 |
| *sfr* | | 输出节点应力破坏比相关的量 |
| *vsf* | | 二维模型时，输出流网 |
| *poreflow* | | 输出孔隙管流计算相关量 |

## DATAPOINT

表 10 DATAPOINT输入格式说明表

|  |  |
| --- | --- |
| **关键词** | **DATAPOINT** |
| 功能 | 输出指定各节点的计算量 |
| 输入格式： | |
| datapoint,*num*=I //命令行  *nnode,issumq,isstat* //控制行  *n1,n2,…,nn //数据行* | |
| 示例： | |
| datapoint,num=1  3 1 0  16 14 12 | |
| 字段说明： | |
| *num* | 输出数据集的个数 |
| *nnode* | 某一数据集的节点数 |
| *issumq* | 是否仅输出数据集各点的流量和  =0，不是，默认  =1，是。 |
| *isstat* | 是否同时输出各量的统计量，包括'SUM', 'MAX', 'X', 'Y', 'Z', 'MIN', 'X', 'Y', 'Z', 'MEAN', 'MEDIAN', 'MAD', 'STD', 'KURTOSIS', 'SKEWNESS'  =0，不是，默认  =1，是。 |
| *n1,n2,…,nn* | 节点号 |

## GE

此参数为基于有限元应力场边坡滑弧进化搜索算法的输入参数。

此命令及RE、WSP,XLIMIT,OPTIMS等6个命令的内容，另存一文件名为“项目名\_eap.dat”的文件，程序会根据项目名自动识别读取。项目名为sinp文件的文件名。

表 11 GE输入格式说明表

|  |  |
| --- | --- |
| **关键词** | **GE** |
| 功能 | 输入边坡坡面线 |
| 输入格式： | |
| GE,*num*=I //命令行  *x1,y1*  *­x2,y2*  *……* | |
| 示例： | |
| GE,NUM=4  0,5  5,5  15,10  25,10 | |
| 字段说明： | |
| *num* | 点数 |
| *xi,yi* | 坡面线点的坐标 |

## RE

此参数为基于有限元应力场边坡滑弧进化搜索算法的输入参数。

此命令及GE、WSP,XLIMIT,OPTIMS等6个命令的内容，另存一文件名为“项目名\_eap.dat”的文件，程序会根据项目名自动识别读取。项目名为sinp文件的文件名。

表 12 RE输入格式说明表

|  |  |
| --- | --- |
| **关键词** | **RE** |
| 功能 | 输入岩面线（底线） |
| 输入格式： | |
| RE,*num*=I //命令行  *x1,y1*  *­x2,y2*  *……* | |
| 示例： | |
| RE,NUM=2  0,0  25,0 | |
| 字段说明： | |
| *num* | 点数 |
| *xi,yi* | 岩面线点的坐标 |

## WSP

此参数为基于有限元应力场边坡滑弧进化搜索算法的输入参数。

此命令及GE、RE, XLIMIT,OPTIMS等6个命令的内容，另存一文件名为“项目名\_eap.dat”的文件，程序会根据项目名自动识别读取。项目名为sinp文件的文件名。

表 12 WSP输入格式说明表

|  |  |
| --- | --- |
| **关键词** | **WSP** |
| 功能 | 输入薄弱夹层底面线 |
| 输入格式： | |
| WSP,*num*=I //命令行  *x1,y1*  *­x2,y2*  *……* | |
| 示例： | |
| WSP,NUM=2  0,12  70,12 | |
| 字段说明： | |
| *num* | 点数 |
| *xi,yi* | 薄弱夹层底面线线点的坐标 |

## XLIMIT

此参数为基于有限元应力场边坡滑弧进化搜索算法的输入参数。

此命令及GE、RE, WSP,OPTIMS等6个命令的内容，另存一文件名为“项目名\_eap.dat”的文件，程序会根据项目名自动识别读取。项目名为sinp文件的文件名。

表 12 XLIMIT输入格式说明表

|  |  |
| --- | --- |
| **关键词** | **XLIMIT** |
| 功能 | 通过输入滑弧进出范围，定义边坡的搜索范围 |
| 输入格式： | |
| xlimit, *xtl*=R, *xtr*=R, *xcl*=R, *xcr*=R | |
| 示例： | |
| xlimit,xtl=0,xtr=5,xcl=15,xcr=25 | |
| 字段说明： | |
| *xtl, xtr* | 滑弧滑出坡脚的最小、最大*x*坐标 |
| *xcl, xcr* | 滑弧滑入坡顶的最小、最大*x*坐标 |

## OPTIMPARA

此参数为基于有限元应力场边坡滑弧进化搜索算法的输入参数。

此命令及GE、RE, WSP,Xlimit等6个命令的内容，另存一文件名为“项目名\_eap.dat”的文件，程序会根据项目名自动识别读取。项目名为sinp文件的文件名。

表 12 2.15 OPTIMPARA输入格式说明表

|  |  |
| --- | --- |
| **关键词** | **OPTIMS/OPTIMPARA** |
| 功能 | 优化算法参数 |
| 输入格式： | |
| optimpara,*Solver=C,nSlice=I,Popsize=I,Maxiter=I,eps1=R,eps2=R,eps3=R,w=R,c1=R,c2=R,gamma=*R*,F=*R*,CR=R,strategy=rand1,sigma=R,mu\_perc=R,Shape=I,IniShape=*I*,nRepeat=I,Icode=I,Imut=I,Irep=I,Icross=I,pcross=R* | |
| 示例： | |
| optimpara,*Solver=GA,nSlice=30,Popsize=60,Maxiter=3000,eps1=0.0001,eps2=0.001,eps3=0.0001,w=0.7298,c1=1.49618,c2=1.49618,gamma=1,F=0.5,CR=0.1,strategy=rand1,sigma=0.5,mu\_perc=0.5,Shape=0,IniShape=1,nRepeat=1,Icode=1,Imut=7,Irep=3,Icross=4,pcross=0.85* | |
| 字段说明： | |
| *solver* | 优化算法，可以为  ga, de(默认), pso, cpso, cmaes |
| *nSlice* | 滑弧分段数  默认30 |
| *popsize* | 初始样本数 |
| *maxiter* | 最大迭代次数 |
| *eps1,eps2,eps3* |  |
| *w,c1,c2,gamma* | PSO算法的参数，  默认值为*w=0.7298,c1=1.49618,c2=1.49618,gamma=1* |
| *F, CR, strategy* | DE算法参数  默认值为*,F=0.5,CR=0.1,strategy=rand1* |
| *sigma, mu\_perc* | cmaes算法参数  默认值为*sigma=0.5,mu\_perc=0.5* |
| *shape* | 滑弧形状  =0，非圆滑弧（默认）  =1，圆弧 |
| *IniShape* | 滑弧的生成算法  =1，圆弧+薄弱层（默认）  =0，cheng氏方法，参考(Cheng YM, Li L, Chi SC. Performance studies on six heuristic global optimization methods in the location of critical slip surface. Computers and Geotechnics. 2007;34(6):462-484.) |
| *nRepeat* | 随机重复计算的次数  =1，默认 |
| *Icode,*  *Imut,*  *Irep, Icross, pcross* | ga算法参数  默认值为*ICODE=1,IMUT=7,ICROSS=4,IREP=3,PCROSS=0.85*   1. icode   0 for integer coded. 1 for real-coded.   1. imut   mutation mode; 1/2/3/4/5 (default for icga is 2).  1=uniform mutation, fixed rate.  2=uniform, adjustable rate based on fitness.  3=uniform adjustable rate based on distance.  4=uniform+creep, fixed rate.  5=uniform+creep, adjustable rate based on fitness.  6=uniform+creep, adjustable rate based on distance.  7=non-uniform mutation for real-coded ga   1. icross   Cross Mode: 0, one-two point x. default for integer-coded'  Cross Mode: 1, for Simple Arithmetic Recombination, for real-coded only'  Cross Mode: 2, for single Arithmetic Recombination, for real-coded only'  Cross Mode: 3, for whole Arithmetic Recombination, for real-coded only'  Cross Mode: 4, for blend-0.5, Recombination, for real-coded only'   1. irep   reproduction plan;  1=Full generational replacement  2=Steady-state-replace-random  3=Steady state-replace-worst (default is 3) elite tournament selection   1. pcross   crossover probability; must be <= 1.0 (default is 0.85).  If crossover takes place, either one or two splicing points are used, with equal probabilities. |

# Exca2D输入文件各关键词输入格式说明

## KPOINT

表 2 KPOINT输入格式说明表

|  |  |  |
| --- | --- | --- |
| **关键词** | **KP** | |
| 功能 | 输入模型关键点信息 | |
| 输入格式： | | |
| Node,*NUM*=I,*DATAPACKING*=0|1,*DIMENSION*=2|3,*isporeflow*=0|1  *X*,*Y*[,*Z*] | | |
| 示例： | | |
| NODE,NUM=16, DATAPACKING=1, DIMENSION=2  0.1465762786596120 0.2506144708076132  0.2500000000000000 0.0000000000000000  …… | | |
| 字段说明： | | |
| *NUM* | | 节点数 |
| *DATAPACKING* | | 输入格式  =1(默认)，按点的顺序输入{x1,y1,z1},{x2,y2,z2},…  =0，按坐标的顺序输入{x1,x2,…},{y1,y2,…},{z1,z3,…} |
| *isporeflow* | | 是否为孔隙管流计算  =0，否（默认）  =1，是 |
| *x,y,z* | | 点的坐标  当dimension=2时，二维模型，输入x,y  当dimension=3时，三维模型，输入x,y,z |