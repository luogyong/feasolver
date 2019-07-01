!���ܣ������ڵ�ˮͷ�仯�������������
!���裬
!1 )�����ʼ��Ӧ������
!2 )����h1���ʼ����ЧӦ����sigma1'��
!3 )�ٶ���Ӧ�����䣬����Dsigma=-(h2-h1)�����ЧӦ������,��sigma2'(=sigma1'+Dsigma')
!4 )ȷ����ж�ؼ�ģ��
!5 ) ����
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

REAL(8) FUNCTION GETU(SELF,IELT,IQU,ISTEP)  !��ȡ��Ԫiel�л��ֵ�IQU,ISTEPʱ��pore water pressure
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
        !�������е��� 
        PH1(1:NNUM1)=PH1(1:NNUM1)+SELF.HEAD(I,ELEMENT(IELT).NODE)*SF(SELF.SF(I)).FACTOR(ISTEP)        
    ENDDO    
    GETU=DOT_PRODUCT(ECP(ET1).LSHAPE(:,IQU),PH1(1:NNUM1))
    GETU=MAX(GETU-ELEMENT(IELT).XYGP(NDIMENSION,IQU),0.0D0)*9.8 !�����ǷǱ��Ͳ����ĸ�ѹ

ENDFUNCTION



endmodule