      !#########################################################################
      !>
      !> \brief GAUSSIAN (NORMAL) CUMULATIVE DISTRIBUTION FUNCTION (USING ERF FINCTION)
      !> \param X ARGUMENT
      !> \param M EXPECTED VALUE
      !> \param S ROOT MEAN SQUARE
      !> \return RETURNS CUMULATIVE DISTRIBUTION FUNCTION
      !>
      REAL FUNCTION GAUSS_CDF(X,M,S)
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      REAL X ! ARGUMENT
      REAL M ! EXPECTED VALUE
      REAL S ! ROOT MEAN SQUARE
      !-------------------------------------------------------------------------
      GAUSS_CDF=(1.+ERF((X-M)/(S*SQRT(2.))))/2.
      RETURN
      END

      !#########################################################################
      !>
      !> \brief GAUSSIAN (NORMAL) CUMULATIVE DISTRIBUTION FUNCTION (USING APPROXIMATION)
      !> \details(APPROXIMATION FROM Research Article High Accurate Simple Approximation
      !> of Normal Distribution Integral. Hector Vazquez-Leal,
      !> Roberto Castaneda-Sheissa, Uriel Filobello-Nino, Arturo Sarmiento-Reyes
      !> Jesus Sanchez Orea)
      !>
      !> \param X ARGUMENT
      !> \param M EXPECTED VALUE
      !> \param S ROOT MEAN SQUARE
      !> \return RETURNS CUMULATIVE DISTRIBUTION FUNCTION
      !>
      REAL FUNCTION GAUSS_CDF_A(X,M,S)
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      REAL X ! ARGUMENT
      REAL M ! EXPECTED VALUE
      REAL S ! ROOT MEAN SQUARE
      !-------------------------------------------------------------------------
      REAL A1,A2,A3
      A1=(X-M)/S
      A2=-358.*A1/23.
      A3=111.*ATAN(37.*A1/294.)
      GAUSS_CDF_A=1./(EXP(A2+A3)+1.)
      RETURN
      END

      !#########################################################################
      !>
      !> \brief CANAL APPROXIMATION OF CHI SQUARED
      !> \details(gives max abs. mistake ~0.008 @ ARG=~0.88, ARG=~2.27)
      !>
      !> \param ARG ARGUMENT
      !> \param V   DEGREE OF FREEDOM
      !> \return RETURNS CDF OF CHI SQUARED
      !>
      REAL FUNCTION CHI_SQ_CANAL(ARG,V)
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      REAL ARG ! ARGUMENT
      REAL V   ! DEGREE OF FREEDOM
      !-------------------------------------------------------------------------
      REAL X,M,S,LOGTV,T
      T=ARG*ARG ! SQUARED !
      LOGTV=LOG(T/V)
      X=EXP((1./6.)*LOGTV)-
     * (EXP((1./3.)*LOGTV)/2.)+
     * (EXP((1./2.)*LOGTV)/3.)
      M=(5./6.)-(1./(9.*V))-(7./(648.*V*V))+(25./(2187.*V*V*V))
      S=SQRT((1./(18.*V))+(3./(162.*V*V))-(37./(11664*V*V*V)))
      !CHI_SQ_CANAL=GAUSS_CDF_A(X,M,S) ! USING APPROXIMATION WITHOUT ERF
      CHI_SQ_CANAL=GAUSS_CDF(X,M,S) ! USING ERF
      RETURN
      END

      !#########################################################################
      !>
      !> \brief KELLEY APPROXIMATION OF CHI SQUARED
      !> \details(gives max abs. mistake ~0.005 @ ARG=~0.4, ARG=~1.27)
      !>
      !> \param ARG ARGUMENT
      !> \param V   DEGREE OF FREEDOM
      !> \return RETURNS CDF OF CHI SQUARED
      !>
      REAL FUNCTION CHI_SQ_KELLEY(ARG,V)
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      REAL ARG ! ARGUMENT
      REAL V   ! DEGREE OF FREEDOM
      !-------------------------------------------------------------------------
      REAL X,T
      T=ARG*ARG ! SQUARED !
      X=(EXP((1./3.)*LOG(T/V))-
     * (1.-(2./(9.*V))) )/
     * SQRT(2./(9.*V))
      !CHI_SQ_KELLEY=GAUSS_CDF_A(X,0.,1.) ! USING APPROXIMATION WITHOUT ERF
      CHI_SQ_KELLEY=GAUSS_CDF(X,0.,1.) ! USING ERF
      RETURN
      END

      !#########################################################################
      !>
      !> \brief PEIZER-PRATT APPROXIMATION OF CHI SQUARED
      !> \details(gives max abs. mistake ~0.0012 @ ARG=~1.74)
      !>
      !> \param ARG ARGUMENT
      !> \param V   DEGREE OF FREEDOM
      !> \return RETURNS CDF OF CHI SQUARED
      !>
      REAL FUNCTION CHI_SQ_PEIZERPRATT(ARG,V)
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      REAL ARG ! ARGUMENT
      REAL V   ! DEGREE OF FREEDOM
      !-------------------------------------------------------------------------
      REAL X,T
      T=ARG*ARG ! SQUARED !
      IF(T.EQ.(V-1))THEN !IF((T-(V-1)).LT.0.001)THEN !
         X=((-1./3.)+(0.08/V))/SQRT((2.*V)-2)
      ELSE
         X=((T-V+(2./3.)-(0.08/V))/ABS(T-(V-1)))*
     *   SQRT((V-1)*LOG((V-1)/T)+T-(V-1))
      ENDIF
      !CHI_SQ_PEIZERPRATT=GAUSS_CDF_A(X,0.,1.) ! USING APPROXIMATION WITHOUT ERF
      CHI_SQ_PEIZERPRATT=GAUSS_CDF(X,0.,1.) ! USING ERF
      RETURN
      END

      !#########################################################################
      !>
      !> \brief POLYNOMIAL APPROXIMATION OF CHI SQUARED WITH 3 DEGREES OF FREEDOM
      !> \details(gives max abs. mistake ~4.5e-5 @ ARG=~5.99)
      !> \author Sobolev A.
      !>
      !> \param ARG ARGUMENT (0..5)
      !> \return RETURNS CDF OF CHI SQUARED
      !>
      REAL FUNCTION CHI_SQ_POLY_3_DEG(ARG)
      !-------------------------------------------------------------------------
      ! LIST OF FORMAL PARAMETERS
      REAL ARG ! ARGUMENT
      !-------------------------------------------------------------------------
      REAL(8) A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,
     *        A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,
     *        A20,A21,A22
      REAL X
      X=ARG
      IF(X.GE.6)THEN
      CHI_SQ_POLY_3_DEG=1.0
      ELSE
      A0 = 3.109446140570875E-012
      A1 = -3.499447493735563E-009
      A2 = 3.986971529314769E-007
      A3 = 2.6595274109083427E-001
      A4 = 8.371102266288947E-005
      A5 = -8.022880310063087E-002
      A6 = 1.4462510174331116E-003
      A7 = 1.105392865463706E-002
      A8 = 4.97525344039133E-003
      A9 = -7.486586198986564E-003
      A10 = 4.7413421552504875E-003
      A11 = -2.7908440930410063E-003
      A12 = 1.392299162863025E-003
      A13 = -4.8671800493258965E-004
      A14 = 1.062998084653591E-004
      A15 = -1.0413499260722634E-005
      A16 = -1.4140096495637238E-006
      A17 = 7.216572747503191E-007
      A18 = -1.3534800869936772E-007
      A19 = 1.5162217924079233E-008
      A20 = -1.0616930719276025E-009
      A21 = 4.323970091673211E-011
      A22 = -7.876584612642333E-013
      CHI_SQ_POLY_3_DEG=REAL
     *    (A0 + X * (A1  + X * (A2  + X * (A3  + X * (A4  + X * (A5
     *        + X * (A6  + X * (A7  + X * (A8  + X * (A9  + X * (A10
     *        + X * (A11 + X * (A12 + X * (A13 + X * (A14 + X * (A15
     *        + X * (A16 + X * (A17 + X * (A18 + X * (A19 + X * (A20
     *        + X * (A21 + X * A22))))))))))))))))))))))
      ENDIF
      RETURN
      END