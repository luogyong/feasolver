MODULE WELL_DS
    
END MODULE
    
    
    
    
REAL(8) Function FrictionDarcy(Model,Epsilon, Diameter, Re, Rew ,PipePoreRatio, RoughFunc, PipeType, Tol) 
    IMPLICIT NONE
    INCLUDE 'DOUBLE.H'
    
  !calculate Darcy friction factor based on Colebrook-White equation
  !Colebrook-White equation
  !1/sqrt(f)=-2.0*lg(epsilon/3.7D+2.51/(Re*sqrt(f)))
  
    CHARACTER(32),INTENT(IN)::Model,PipeType
    REAL(DPN),INTENT(IN)::Epsilon,Diameter,Re,Rew, PipePoreRatio, RoughFunc,Tol
    REAL(DPN)::X,F,DF,DX,FOU,F0,FT,Dke,Epsilon1,T1
    INTEGER::I


 
    Dke=0
    If(Model=="Siwon".And.PipePoreRatio>0.01)Dke=Diameter*0.282*PipePoreRatio**2.4
    Epsilon1=Epsilon+Dke

    If(Re<=2320)Then

        FrictionDarcy=64.0/Re

    Else
        !approximatedwithSwamee-Jainequation
        x=2.0*Abs(Log10(Epsilon1/(3.7*Diameter)+5.74/Re**0.9))
        dx=0.
        i=0
        Do While((Abs(dx)>tol).And.(i<101))
            i=i+1
            x=x+dx
            f=x+2.0*(Log10(Epsilon1/(3.7*Diameter)+2.51*x/Re))
            df=1+2.1801583/(Epsilon1*Re/(3.7*Diameter)+2.51*x)
            dx=-f/df
        ENDDO
        FrictionDarcy=x**(-2.0)
        !'At(0)=At(0)*(1+0.0835)'justforSuTest1dataonly
        !At(1)=i
    EndIf


    If(i>100) STOP "Failedincalculatingdarcy(normal)friction."
        
    

    Select Case(Model)

    Case("Ouyang")
        !MsgBox("ModelisOuyang.")
        If (Re<=2320) Then
            If(Rew>=0.0)Then
                fou=(1+0.04304*Rew**0.6142)
            Else
                fou=(1-0.0625*(-Rew)**1.3056/(Rew+4.626)**(-0.2724))
            EndIf
        Else
            If(Rew>=0.0)Then
                If (PipeType=="Wellbore") fou=(1-0.0153*Rew**0.3978) !wellbore
                If(PipeType=="Porousbore") fou=(1-29.03*(Rew/Re)**0.8003) !'porouspipe
            Else
                fou=(1-17.5*Rew/Re**0.75)
            EndIf
        EndIf

        FrictionDarcy=FrictionDarcy*fou

    Case("Siwon")
        f0=0.0106*PipePoreRatio**0.413
        FrictionDarcy=FrictionDarcy+f0
    Case("SuZe","WangZM")

        x=FrictionDarcy

        If(Abs(RoughFunc)>0.000001)Then
            t1=(8.0/x)**0.5-1.25*Log(x)-RoughFunc
            dx=0.0
            i=0
            Do While ((Abs(dx)>tol).And.(i<101))
                i=i+1
                x=x+dx
                f=(8.0/x)**0.5-1.25*Log(x)-t1
                df=-0.5*8**0.5*x**(-1.5)-1.25/x
                dx=-f/df
            end do
            If(i>100)  STOP "FailedincalculateSuZeorWangZMfriction."
        EndIf



        FrictionDarcy=x
        !!!!!!!'''''
        If(Model=="WangZM")  FrictionDarcy=x*(1-PipePoreRatio)


    Case("ChenYS")

    Case DEFAULT
    !i=MsgBox("NOSUCHMODELNAMED"&Model&Chr(13)&"ThenormalfDisReturned.",0)
    !Ifi=1ThenExitFunction

    EndSelect

End Function

REAL(8) Function Factor_ACC_HeadLoss_Term(Model,PipePoreRatio,VraRatio)
    !VraRatio= inflow velocity / main flow velocity =Vw/Va=phi*Vbp/Vb
    IMPLICIT NONE
    INCLUDE 'DOUBLE.H'
    
    CHARACTER(32),INTENT(IN)::Model
    REAL(DPN),INTENT(IN)::PipePoreRatio,VraRatio 
    REAL(DPN)::FC1,B1,BETA1,T1
    
    Select Case(Model)
    
    Case("Siwon")

        If(VraRatio<0.002)Then
            beta1=1.024
        ElseIf(VraRatio>0.017)Then
            beta1=1.11
        Else
            beta1=1.034+4.27*VraRatio
        EndIf

        !beta1=1.05

        If(Abs(VraRatio)<1.0D-10)Then
            fc1=beta1
        Else
            b1=10.0/(1000*PipePoreRatio)**4.2+4.0/10**7
            T1=VraRatio/PipePoreRatio !Vbp/Vb
            fc1=beta1*(1+1.235/beta1/(b1/(T1)**2+1.235)**2.0)
        EndIf

    Case("WangZM")

        If(VraRatio<0.002)Then
            beta1=1.024
        ElseIf(VraRatio>0.017)Then
            beta1=1.11
        Else
            beta1=1.034+4.27*VraRatio
        EndIf

        !beta1=1.0+5.0*VraRatio

        fc1=2.0*beta1

    Case DEFAULT
            fc1=2.0
    EndSelect

    Factor_ACC_HeadLoss_Term=fc1
    

End Function

REAL(8) Function RoughnessFunction(Model,Re,PipePoreRatio,VraRatio,PerforationDensity,dDRatio)
    IMPLICIT NONE
    INCLUDE 'DOUBLE.H'
    CHARACTER(32),INTENT(IN)::Model
    REAL(DPN),INTENT(IN)::Re,PipePoreRatio,VraRatio,PerforationDensity,dDRatio
    
    If (Model == "WangZM") Then
        RoughnessFunction = Log(Re) * (PipePoreRatio - VraRatio)
    
    ElseIf (Model == "SuZe") Then
        RoughnessFunction = 0.1778 * PerforationDensity * dDRatio
    Else
        RoughnessFunction = 0.0
    End If

End Function

REAL(8) Function WaterDensity(T)
    IMPLICIT NONE
    INCLUDE 'DOUBLE.H'
    REAL(DPN),INTENT(IN)::T
    REAL(DPN)::C(10)
    
    C(1)=999.83952
    C(2)=16.945176
    C(3)=-0.0079870401
    C(4)=-0.000046170461
    C(5)=0.00000010556302
    C(6)=-2.8054253E-10
    C(7)=1.2378
    C(8)=-0.001303
    C(9)=0.00000306
    C(10)=0.0000000255

    WaterDensity=(C(1)+C(2)*T+C(3)*T**2+C(4)*T**3+C(5)*T**4+C(6)*T**5)/(1+0.01687985*T)
End Function

REAL(8) Function DynamicVisCosity(DV20,T)
    
    IMPLICIT NONE
    INCLUDE 'DOUBLE.H'
    
    REAL(DPN),INTENT(IN)::DV20,T
    REAL(DPN)::C(10)
    
    C(1)=999.83952
    C(2)=16.945176
    C(3)=-0.0079870401
    C(4)=-0.000046170461
    C(5)=0.00000010556302
    C(6)=-2.8054253E-10
    C(7)=1.2378
    C(8)=-0.001303
    C(9)=0.00000306
    C(10)=0.0000000255

    DynamicVisCosity=(20.0-T)/(96.0+T)*(C(7)+C(8)*(20.0-T)+C(9)*(20.0-T)**2.0+C(10)*(20.0-T)**3.0)
    DynamicVisCosity=DV20*10.0**DynamicVisCosity
    
    End Function
    
    