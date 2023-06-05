
Modifications comparied with the CosmoMC-planck2018 package.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
file camb/equations_ppfi.f90
	add module IDEtools

    add module CoupledFluidModels

    add module CoupledQuintModels

    rewrite module LambdaGeneral
    
    rewrite subroutine init_background

    modify function dtauda(a)
        ..................................................................
        integer nu_i
        real(dl) grhodea4,grhoca4 !8*pi*G*rho_{de}*a**4, 8*pi*G*rho_{c}*a**4, add for IDE
        a2=a**2

        !  8*pi*G*rho*a**4.
        !grhoa2=grhok*a2+(grhoc+grhob)*a+grhog+grhornomass
        if(a==0._dl) then
            grhodea4=0._dl
            grhoca4=0._dl
        else
            call IDEout(a,grhodea4,grhoca4)
            grhodea4=grhodea4*a2
            grhoca4=grhoca4*a2
        end if 
        grhoa2=grhok*a2+grhob*a+grhog+grhornomass+grhoca4+grhodea4 ! Modified due to dark sector interaction

        if (CP%Num_Nu_massive /= 0) then
        ..................................................................

    modify subroutine SetupScalarArrayIndices(EV, max_num_eqns)
        ..................................................................
        maxeq = maxeq +  (EV%lmaxg+1)+(EV%lmaxnr+1)+EV%lmaxgpol-1

        ! !Dark energy
        ! if (.not. is_cosmological_constant) then
            ! EV%w_ix = neq+1
            ! neq=neq+1 !ppf
            ! maxeq=maxeq+1
        ! else
            ! EV%w_ix=0
        ! end if
        if (w_Perturb) then
            if (use_ppf) then
                EV%w_ix = neq+1
                neq=neq+1 !ppf
                maxeq=maxeq+1
                if (evolve_vc) then
                   neq=neq+1 !ide vc
                   maxeq=maxeq+1
                end if
            else
                EV%w_ix = neq+1
                neq=neq+2
                maxeq=maxeq+2
                if (evolve_vc) then
                   neq=neq+1 !ide vc
                   maxeq=maxeq+1
                end if
            end if
        else
            EV%w_ix=0
        end if

        !Massive neutrinos
        ..................................................................

    modify subroutine CopyScalarVariableArray(y,yout, EV, EVout)
        ..................................................................
        yout(1:basic_num_eqns) = y(1:basic_num_eqns)
        ! if (.not. is_cosmological_constant) then
            ! yout(EVout%w_ix)=y(EV%w_ix)
        ! end if
        if (w_Perturb) then
            if (use_ppf) then
                yout(EVout%w_ix)=y(EV%w_ix)
                if (evolve_vc) then  !add for vc 
                    yout(EVout%w_ix+1)=y(EV%w_ix+1)
                end if
            else
                yout(EVout%w_ix)=y(EV%w_ix)
                yout(EVout%w_ix+1)=y(EV%w_ix+1)
                if (evolve_vc) then  !add for vc 
                    yout(EVout%w_ix+2)=y(EV%w_ix+2)
                end if
            end if
        end if

        if (.not. EV%no_phot_multpoles .and. .not. EVout%no_phot_multpoles) then
        ..................................................................

    modify subroutine output(EV,y, tau,sources, notused_custom_sources)
        ..................................................................
        real(dl) opacity, dopacity, ddopacity, visibility, dvisibility, ddvisibility, exptau, lenswindow
        real(dl) vc, clxq, vq !for IDE

        call IonizationFunctionsAtTime(tau, opacity, dopacity, ddopacity, &
        visibility, dvisibility, ddvisibility, exptau, lenswindow)
        ..................................................................
        vbdot =yprime(5)

        !  Compute expansion rate from: grho 8*pi*rho*a**2

        grhob_t=grhob/a
        !grhoc_t=grhoc/a
        call IDEout(a,grhov_t,grhoc_t,wde=w_eff) 
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
    
        !  8*pi*a*a*SUM[rho_i*clx_i] add radiation later
        dgrho=grhob_t*clxb+grhoc_t*clxc
    
        !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
        !dgq=grhob_t*vb
        !add for IDE
        vc=0
        if (w_Perturb .and. evolve_vc) then
            if (use_ppf) then
                vc=y(EV%w_ix+1)
            else
                vc=y(EV%w_ix+2)
            end if
        end if
        dgq=grhob_t*vb+grhoc_t*vc
    
        ! if (is_cosmological_constant) then
            ! w_eff = -1_dl
            ! grhov_t=grhov*a2
        ! else
            ! !ppf
            ! w_eff=w_de(a)   !effective de
            ! grhov_t=grho_de(a)/a2
            ! dgrho=dgrho+EV%dgrho_e_ppf
            ! dgq=dgq+EV%dgq_e_ppf
        ! end if
        !add for IDE
        if (w_Perturb) then
            if (use_ppf) then
                dgrho=dgrho+EV%dgrho_e_ppf
                dgq=dgq+EV%dgq_e_ppf
            else
                clxq=y(EV%w_ix)
                vq=y(EV%w_ix+1)
                dgrho=dgrho + clxq*grhov_t
                dgq = dgq + vq*grhov_t*(1+w_eff)
            end if
        end if
        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        ..................................................................
        sigma=(z+1.5_dl*dgq/k2)/EV%Kf(1)

        ! if (is_cosmological_constant) then
            ! ppiedot=0
        ! else
            ! hdotoh=(-3._dl*grho-3._dl*gpres -2._dl*grhok)/6._dl/adotoa
            ! ppiedot=3._dl*EV%dgrho_e_ppf+EV%dgq_e_ppf*(12._dl/k*adotoa+k/adotoa-3._dl/k*(adotoa+hdotoh))+ &
                ! grhov_t*(1+w_eff)*k*z/adotoa -2._dl*k2*EV%Kf(1)*(yprime(EV%w_ix)/adotoa-2._dl*y(EV%w_ix))
            ! ppiedot=ppiedot*adotoa/EV%Kf(1)
        ! end if
        ppiedot=0
        if (w_Perturb .and. use_ppf) then
            hdotoh=(-3._dl*grho-3._dl*gpres -2._dl*grhok)/6._dl/adotoa
            ppiedot=3._dl*EV%dgrho_e_ppf+EV%dgq_e_ppf*(12._dl/k*adotoa+k/adotoa-3._dl/k*(adotoa+hdotoh))+ &
                grhov_t*(1+w_eff)*k*z/adotoa -2._dl*k2*EV%Kf(1)*(yprime(EV%w_ix)/adotoa-2._dl*y(EV%w_ix))
            ppiedot=ppiedot*adotoa/EV%Kf(1)
        end if

        polter = 0.1_dl*pig+9._dl/15._dl*ypol(2)
        ..................................................................

    modify subroutine initial(EV,y, tau)
        ..................................................................
        y(EV%g_ix+1)=InitVec(i_qg)

        ! if (.not. is_cosmological_constant) then
            ! y(EV%w_ix) = InitVec(i_clxq) !ppf: Gamma=0, i_clxq stands for i_Gamma
        ! end if
        if (w_Perturb) then
            if (use_ppf) then
                y(EV%w_ix) = InitVec(i_clxq) !ppf: Gamma=0, i_clxq stands for i_Gamma
                if (evolve_vc) then
                    y(EV%w_ix+1) = InitVec(i_vq) !ide: vc=0, i_vq stands for i_vc
                end if
            else 
                y(EV%w_ix) = InitVec(i_clxq)
                y(EV%w_ix+1) = InitVec(i_vq)
                if (evolve_vc) then
                    y(EV%w_ix+2) = InitVec(i_vq) !ide: vc=0, i_vq stands for i_vc
                end if
            end if
        end if

        !  Neutrinos
        y(EV%r_ix)=InitVec(i_clxr)
        ..................................................................

    modify subroutine derivs(EV,n,tau,ay,ayprime) 
        ..................................................................
        real(dl) w_eff, grhoT       
        !IDE
        real(dl) grhodea4,grhoca4,gQ,gC(3),gD(2),gC1,gC2,gC3,gD1,gD2,S0,S1,S2,xi0,deltapT,vc,vcdot,ca2

        k=EV%k_buf
        ..................................................................
        !  Compute expansion rate from: grho 8*pi*rho*a**2

        grhob_t=grhob/a
        !grhoc_t=grhoc/a
        call IDEout(a,grhov_t,grhoc_t,wde=w_eff,gQ=gQ,gC=gC,gD=gD,ca2=ca2)
        gC1=gC(1); gC2=gC(2); gC3=gC(3); gD1=gD(1); gD2=gD(2) 
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        ! if (is_cosmological_constant) then
            ! grhov_t=grhov*a2
            ! w_eff = -1_dl
        ! else
            ! !ppf
            ! w_eff=w_de(a)   !effective de
            ! grhov_t=grho_de(a)/a2
        ! end if

        !  Get sound speed and ionisation fraction.
        ..................................................................
        !total perturbations: matter terms first, then add massive nu, de and radiation
        !  8*pi*a*a*SUM[rho_i*clx_i]
        dgrho_matter=grhob_t*clxb+grhoc_t*clxc
        !  8*pi*a*a*SUM[(rho_i+p_i)*v_i]
        !dgq=grhob_t*vb
        vc=0
        if (w_Perturb .and. evolve_vc) then
            if (use_ppf) then
                vc=ay(EV%w_ix+1)
            else
                vc=ay(EV%w_ix+2)
            end if
        end if	
        dgq=grhob_t*vb + grhoc_t*vc

        if (CP%Num_Nu_Massive > 0) then
        ..................................................................
        ! if (w_lam /= -1 .and. w_Perturb) then
        !    clxq=ay(EV%w_ix)
        !    vq=ay(EV%w_ix+1)
        !    dgrho=dgrho + clxq*grhov_t
        !    dgq = dgq + vq*grhov_t*(1+w_lam)
        !end if
        if (w_Perturb .and. .not. use_ppf) then
            clxq=ay(EV%w_ix)
            vq=ay(EV%w_ix+1)
            dgrho=dgrho + clxq*grhov_t
            dgq = dgq + vq*grhov_t*(1+w_eff)
        end if

        if (EV%no_nu_multpoles) then
        ..................................................................
        ! if (.not. is_cosmological_constant) then
            ! !ppf
            ! grhoT = grho - grhov_t
            ! vT= dgq/(grhoT+gpres)
            ! Gamma=ay(EV%w_ix)

            ! !sigma for ppf
            ! sigma = (etak + (dgrho + 3*adotoa/k*dgq)/2._dl/k)/EV%kf(1) - k*Gamma
            ! sigma = sigma/adotoa

            ! S_Gamma=grhov_t*(1+w_eff)*(vT+sigma)*k/adotoa/2._dl/k2
            ! ckH=c_Gamma_ppf*k/adotoa
            ! Gammadot=S_Gamma/(1+ckH*ckH)- Gamma -ckH*ckH*Gamma
            ! Gammadot=Gammadot*adotoa
            ! ayprime(EV%w_ix)=Gammadot

            ! if(ckH*ckH.gt.3.d1)then
                ! Gamma=0
                ! Gammadot=0.d0
                ! ayprime(EV%w_ix)=Gammadot
            ! endif

            ! Fa=1+3*(grhoT+gpres)/2._dl/k2/EV%kf(1)
            ! dgqe=S_Gamma - Gammadot/adotoa - Gamma
            ! dgqe=-dgqe/Fa*2._dl*k*adotoa + vT*grhov_t*(1+w_eff)
            ! dgrhoe=-2*k2*EV%kf(1)*Gamma-3/k*adotoa*dgqe
            ! dgrho=dgrho+dgrhoe
            ! dgq=dgq+dgqe

            ! EV%dgrho_e_ppf=dgrhoe
            ! EV%dgq_e_ppf=dgqe
        ! end if
        if (w_Perturb .and. use_ppf) then
            grhoT = grho - grhov_t
            vT= dgq/(grhoT+gpres)
            Gamma=ay(EV%w_ix)
            if(abs(1._dl+w_eff)<1.d-8) then
                ckH=c_Gamma_ppf*k/adotoa
                Gammadot = -(ckH*ckH+2)*adotoa*Gamma
                ayprime(EV%w_ix)=Gammadot
                if (ckH*ckH.gt.3.d1) then
                    Gamma=0
                    Gammadot=0.d0
                    ayprime(EV%w_ix)=Gammadot
                endif
                dgqe=0._dl
                dgrhoe=-2*k2*EV%kf(1)*Gamma + gQ*vT/k
                dgrho=dgrho+dgrhoe
                dgq=dgq+dgqe
                clxq=dgrhoe/grhov_t
                vq=0._dl
                EV%dgrho_e_ppf=dgrhoe
                EV%dgq_e_ppf=dgqe
            else
            !sigma for ppf
            !sigma = (etak + (dgrho + 3*adotoa/k*dgq)/2._dl/k)/EV%kf(1) - k*Gamma
            !Modified due to dark sector interaction 
            if(grhov_t==0._dl .or. abs(a-a_trans)/a<0.02 ) then
                dgqe=0._dl
                dgrhoe=0._dl
                dgrho=dgrho+dgrhoe
                dgq=dgq+dgqe

                clxq=0._dl
                vq=0._dl

                EV%dgrho_e_ppf=dgrhoe
                EV%dgq_e_ppf=dgqe
            else
            sigma = (etak + (dgrho + 3*adotoa/k*dgq + gQ*vT/k)/2._dl/k)/EV%kf(1) - k*Gamma 
            sigma = sigma/adotoa
            dgpi = grhor_t*pir + grhog_t*pig
            if (CP%Num_Nu_Massive > 0) then
                call MassiveNuVarsOut(EV,ay,ayprime,a,dgpi)
            end if
            deltapT=(grhog_t*clxg+grhor_t*clxr+4*(grhog_t+grhor_t)*vT*adotoa/k)/3._dl
            xi0 = -(deltapT - 2*EV%kf(1)*dgpi/3._dl - gD2*(vc-vT)/k)/(grhoT+gpres)
            Fa=1+3*(grhoT+gpres)/2._dl/k2/EV%kf(1)

            S0 = -3*gD2*(vc-vT)*adotoa/k - gC2*(clxc +(3*adotoa+gQ/grhoc_t)*vT/k) - gQ*xi0
            S0 = grhov_t*(1+w_eff)*(vT+sigma)*k + S0/EV%kf(1)  
            S0 = S0/2._dl/k2
            ckH=c_Gamma_ppf*k/adotoa
            Gammadot=S0/(ckH*ckH+1)-(adotoa-gC1/grhov_t)*(ckH*ckH+1)*Gamma
            ayprime(EV%w_ix)=Gammadot
            if (ckH*ckH.gt.3.d1) then
                Gamma=0
                Gammadot=0.d0
                ayprime(EV%w_ix)=Gammadot
            endif
            S_Gamma = S0+gC1/grhov_t*Gamma
            dgqe=S_Gamma - Gammadot - Gamma*adotoa
            dgqe=-dgqe/Fa*2._dl*k + vT*grhov_t*(1+w_eff)
            dgrhoe=-2*k2*EV%kf(1)*Gamma-3/k*adotoa*dgqe + gQ*vT/k

            dgrho=dgrho+dgrhoe
            dgq=dgq+dgqe

            clxq=dgrhoe/grhov_t
            vq=dgqe/grhov_t/(1+w_eff)


            EV%dgrho_e_ppf=dgrhoe
            EV%dgq_e_ppf=dgqe
        end if
        end if
        end if

        !  Get sigma (shear) and z from the constraints
        ..................................................................
        !if (w_lam /= -1 .and. w_Perturb) then
        !
        !   ayprime(EV%w_ix)= -3*adotoa*(cs2_lam-w_lam)*(clxq+3*adotoa*(1+w_lam)*vq/k) &
        !       -(1+w_lam)*k*vq -(1+w_lam)*k*z
        !
        !   ayprime(EV%w_ix+1) = -adotoa*(1-3*cs2_lam)*vq + k*cs2_lam*clxq/(1+w_lam)
        !
        !end if
        !
        if (w_Perturb .and. .not. use_ppf) then
            ayprime(EV%w_ix)= -3*adotoa*(cs2_lam-w_eff)*clxq-9*adotoa**2*(cs2_lam-ca2)*(1+w_eff)*vq/k-(1+w_eff)*k*(vq+z) &
                +gQ/grhov_t*(-clxq+3*adotoa*(cs2_lam-ca2)*vq/k)+(gC1*clxq+gC2*clxc+gC3*vq/k)/grhov_t

            ayprime(EV%w_ix+1) = -adotoa*(1-3*cs2_lam)*vq + k*cs2_lam*clxq/(1+w_eff) &
                    +(gD1*vq+gD2*vc-gQ*(1+cs2_lam)*vq)/grhov_t/(1+w_eff)
        end if

        if (associated(EV%OutputTransfer)) then
        ..................................................................
        !  CDM equation of motion
        !clxcdot=-k*z
        clxcdot=-k*z-k*vc+(gQ*clxc-gC1*clxq-gC2*clxc-gC3*vq/k)/grhoc_t
        ayprime(3)=clxcdot
        if (w_Perturb .and. evolve_vc) then
            vcdot = -adotoa*vc+(gQ-gD2)/grhoc_t*vc-gD1/grhoc_t*vq
            if (use_ppf) then
                ayprime(EV%w_ix+1) = vcdot
            else
                ayprime(EV%w_ix+2) = vcdot
            end if
        end if

        !  Baryon equation of motion.
        clxbdot=-k*(z+vb)
        ..................................................................

    modify subroutine derivsv(EV,n,tau,yv,yvprime)
        ..................................................................
        real(dl) pir,adotoa
        !for IDE
        real(dl) w_eff

        stop 'ppf not implemented for vectors'
        ..................................................................
        ! Compute expansion rate from: grho=8*pi*rho*a**2
        ! Also calculate gpres: 8*pi*p*a**2
        grhob_t=grhob/a
        !grhoc_t=grhoc/a
        call IDEout(a,grhov_t,grhoc_t,wde=w_eff)
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        !grhov_t=grhov*a**(-1-3*w_lam)

        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        !gpres=(grhog_t+grhor_t)/3._dl+grhov_t*w_lam
        gpres=(grhog_t+grhor_t)/3._dl+grhov_t*w_eff
        ..................................................................

    modify subroutine derivst(EV,n,tau,ayt,aytprime)
        ..................................................................
        ! Compute expansion rate from: grho=8*pi*rho*a**2
        ! Also calculate gpres: 8*pi*p*a**2
        grhob_t=grhob/a
        !grhoc_t=grhoc/a
        call IDEout(a,grhov_t,grhoc_t) 
        grhor_t=grhornomass/a2
        grhog_t=grhog/a2
        ! if (is_cosmological_constant) then
            ! grhov_t=grhov*a2
        ! else
            ! grhov_t=grho_de(a)/a2
        ! end if

        grho=grhob_t+grhoc_t+grhor_t+grhog_t+grhov_t
        ..................................................................

file source/CosmologyTypes.f90
    modify Type TCosmoTheoryParams
	    ..................................................................
        integer :: neutrino_hierarchy = neutrino_hierarchy_normal	
		
		!for IDE
        logical :: Run_H0 = .false.
        logical :: Use_PPF = .true.
        integer :: NP_IDE = 0
        integer :: Class_IDE = -1 !any value less than 1 corresponding to default non-coupled CPL model 
        integer :: WForm_CF = 1, QForm_CF = 1, CovQForm_CF = 1
        integer :: UForm_CQ = 1, QForm_CQ = 1 
        character(LEN=1024) :: paramnames_filename = 'paramnames/params_CMB.paramnames'
		character(LEN=1024) :: test_output_root = ''
        ..................................................................

    modify Type, extends(TTheoryParams) :: CMBParams
        ..................................................................
        real(mcp) w, wa
        real(mcp) w0, w1, c_hde, beta_cf, alpha_quint, beta_cq  !ide parameters
        real(mcp) YHe, nnu, iso_cdm_correlated, ALens, Alensf, fdm !fdm is dark matter annihilation, eg,. 0910.3663
        ..................................................................

    modify subroutine TCosmoTheorySettings_ReadParams(this, Ini)
        ..................................................................
        call Ini%Read('lmax_tensor',this%lmax_tensor)

        !for IDE
        this%Run_H0 = Ini%Read_Logical('Run_H0',.false.)
        this%Class_IDE = Ini%Read_Int('Class_IDE',-1)
        
        if (Ini%HasKey('Use_PPF')) then
            call Ini%Read('Use_PPF',this%Use_PPF)
            write(*,*) 'Use_PPF now set from ini file, and its value is: ', this%Use_PPF
        else
            this%Use_PPF=.true.
            if(this%Class_IDE==2) this%Use_PPF=.false. !if DE is quintessence, not use ppf approach.
            write(*,*) 'NOTE: Use_PPF set internally, whose value depends on the specific IDE Class: ',this%Use_PPF
        end if
        
        select case(this%Class_IDE)
        case (1)
            this%WForm_CF = Ini%Read_Int('WForm_CF',1)
            this%QForm_CF = Ini%Read_Int('QForm_CF',1)
            this%CovQForm_CF = Ini%Read_Int('CovQForm_CF',1)
            this%NP_IDE = Ini%Read_Int('NP_CF',0)
            this%paramnames_filename = Ini%ReadFilename('paramnames_CF')
            write(*,*) 'The cosmological model is coupled fluid model'
            write(*,*) 'WForm, QForm, CovQForm is: ', this%WForm_CF,this%QForm_CF,this%CovQForm_CF
        case (2)
            this%UForm_CQ = Ini%Read_Int('UForm_CQ',1)
            this%QForm_CQ = Ini%Read_Int('QForm_CQ',1)
            this%NP_IDE = Ini%Read_Int('NP_CQ',0)
            this%paramnames_filename = Ini%ReadFilename('paramnames_CQ')
            write(*,*) 'The cosmological model is coupled quintessence model'
            write(*,*) 'UForm, QForm is: ', this%UForm_CQ,this%QForm_CQ
        case default
            this%NP_IDE = 0
            this%paramnames_filename = 'paramnames/params_CMB.paramnames'
            write(*,*) 'The cosmological model is the default model'
        end select
        if(Ini%Read_Int('action')==4) this%test_output_root = Ini%Read_String('test_output_root')
        ..................................................................

file source/Calculator_CAMB.f90

    modify subroutine CAMBCalc_CMBToCAMB(this,CMB,P)
        ..................................................................
        wa_ppf = CMB%wa

        !add parameter for IDE
        select case(CosmoSettings%Class_IDE)
        case (1)
            CFP%w0 = CMB%w0
            CFP%w1 = CMB%w1
            CFP%c = CMB%c_hde
            CFP%beta = CMB%beta_cf
        case (2)
            CQP%alpha = CMB%alpha_quint
            CQP%beta = CMB%beta_cq
        end select

        ALens = CMB%ALens
        ..................................................................

    modify subroutine CAMBCalc_InitCAMBParams(this,P)
        ..................................................................
        use ModelParams
        use LambdaGeneral, only : use_ppf, IDE_Class,CFT,CQT,test_output_root
        class(CAMB_Calculator) :: this
        ..................................................................
        P%Transfer%k_per_logint=0

        !for ide
        IDE_Class = CosmoSettings%Class_IDE
        use_ppf = CosmoSettings%Use_PPF
        test_output_root = CosmoSettings%test_output_root
        if (CosmoSettings%Class_IDE>0) then     
            CFT%WForm = CosmoSettings%WForm_CF
            CFT%QForm = CosmoSettings%QForm_CF 
            CFT%CovQForm = CosmoSettings%CovQForm_CF 
            CQT%UForm = CosmoSettings%UForm_CQ
            CQT%QForm = CosmoSettings%QForm_CQ 
        end if

        if (CosmoSettings%use_nonlinear) then
        ..................................................................

file source/CosmologyParameterizations.f90
    modify subroutine TP_Init(this, Ini, Names, Config)
        ..................................................................
        !call this%Initialize(Ini,Names, 'paramnames/params_CMB.paramnames', Config)
        call this%Initialize(Ini,Names, CosmoSettings%paramnames_filename, Config) !for IDE
        if (CosmoSettings%bbn_consistency) call Names%Add('paramnames/derived_bbn.paramnames')
        ..................................................................
        !set number of hard parameters, number of initial power spectrum parameters
        !call this%SetTheoryParameterNumbers(16,last_power_index)
        call this%SetTheoryParameterNumbers(16+CosmoSettings%NP_IDE,last_power_index) !changed for the ide models
        ..................................................................

    modify subroutine TP_ParamArrayToTheoryParams(this, Params, CMB)
        ..................................................................
        call this%TCosmologyParameterization%ParamArrayToTheoryParams(Params, CMB)

        !directly use H0 as free param for IDE
        if(CosmoSettings%Run_H0) then
            DA = Params(3) !directly use H0 as free param
            call SetForH(Params,CMB,DA, .true.)
            !!call InitCAMB(CMB,error)
            if (CMB%tau==0._mcp) then
                CMB%zre=0
            else
                CMB%zre = CosmoCalc%GetZreFromTau(CMB, CMB%tau)
            end if

            LastCMB(cache) = CMB
            cache = mod(cache,ncache)+1
            return
        end if
        
        error = 0   !JD to prevent stops when using bbn_consistency or m_sterile
        ..................................................................

    modify subroutine SetForH(Params,CMB,H0, firsttime,error)
        ..................................................................
        CMB%fdm = Params(16)

        !for IDE
        select case(CosmoSettings%Class_IDE)
        case (1)
            CMB%w0 = Params(17)
            CMB%w1 = Params(18)
            CMB%c_hde = Params(19)
            CMB%beta_cf = Params(20)
        case (2)
            CMB%alpha_quint = Params(17)
            CMB%beta_cq = Params(18)
        end select
        
        call SetFast(Params,CMB)
        ..................................................................

file paramnames/params_CF.paramnames added

file paramnames/params_CQ.paramnames added

file test_ide.ini added

file source/Makefile modified

	#RECOMBINATION=$(RECOMBINATION) EQUATIONS=equations_ppf NONLINEAR=halofit_ppf
	RECOMBINATION=$(RECOMBINATION) EQUATIONS=equations_ppfi NONLINEAR=halofit_ppf
	


