!功能，给出节点水头变化量，算出沉降。
!步骤，
!1 )算出初始总应力场；
!2 )根据h1算初始的有效应力场sigma1'；
!3 )假定总应力不变，根据Dsigma=-(h2-h1)算出有效应力增量,及sigma2'(=sigma1'+Dsigma')
!4 )确定加卸载及模量
!5 ) 计算
module DownWaternSettlement
USE solverds,ONLY:NDIMENSION,NODE,NNUM,ELEMENT,ENUM,ECP,SF
implicit none

public::settlement_head

private

type settlement_load_tydef
    integer::nhead=0,nnode=0
    real(8),allocatable::head(:,:)
    integer,allocatable::sf(:)
contains
    procedure::getpwp=>getu
endtype
type(settlement_load_tydef)::settlement_head

contains

REAL(8) FUNCTION GETU(SELF,IELT,IQU,ISTEP)  !获取单元iel中积分点IQU,ISTEP时的pore water pressure
    CLASS(settlement_load_tydef)::SELF
    INTEGER,INTENT(IN)::IELT,IQU,ISTEP !itime=0,initial, others,changed
    INTEGER::I,ET1,NNUM1,it1
    REAL(8)::PH1(50)
    
    GETU=0.D0
    IF(SELF.NNODE<=0) RETURN
    
    NNUM1=ELEMENT(IELT).NNUM
    ET1=ELEMENT(IELT).ET
    PH1(1:NNUM1)=0.D0
    DO I=1,SELF.NHEAD
        IF(ABS(SF(SELF.SF(I)).FACTOR(ISTEP)+999)<1.E-7) CYCLE
        !各步进行叠加 
        PH1(1:NNUM1)=PH1(1:NNUM1)+SELF.HEAD(I,ELEMENT(IELT).NODE)*SF(SELF.SF(I)).FACTOR(ISTEP)        
    ENDDO    
    GETU=DOT_PRODUCT(ECP(ET1).LSHAPE(:,IQU),PH1(1:NNUM1))
    GETU=MAX(GETU-ELEMENT(IELT).XYGP(NDIMENSION,IQU),0.0D0)*9.8 !不考虑非饱和产生的负压

ENDFUNCTION



endmodule