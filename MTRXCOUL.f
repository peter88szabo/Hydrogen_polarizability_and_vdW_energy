CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      FUNCTION RINIT(N,L,X,Z)
C
C   RADIAL HYDROGENIC WAVE FUNCTIONS R_NL(X,Z) FOR BOUND STATES.
C   N: PRINCIPAL QUANTUM NUMBER
C   L: ANGULAR MOMENTUM QUANTUM NUMBER
C   X=Z*R, WHERE:
C        Z: NUCLEAR CHARGE
C        R: RADIAL DISTANCE IN ATOMIC UNITS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      ANORM=DSQRT(Z)*Z
      RINIT=0.D0
      IF(X.GT.150.D0)RETURN
      GOTO(1,2,3)N
1     RINIT=2.D0*ANORM*DEXP(-X)
      RETURN
2     IF(L.EQ.0)RINIT=ANORM/DSQRT(8.D0)*(2.D0-X)*DEXP(-.5D0*X)
      IF(L.EQ.1)RINIT=ANORM/DSQRT(24.D0)*X*DEXP(-.5D0*X)
      RETURN
3     IF(L.EQ.0)RINIT=ANORM/9.D0/DSQRT(3.D0)*
     * (6.D0-4.D0*X+4.D0/9.D0*X**2)*DEXP(-X/3.D0)
      IF(L.EQ.1)RINIT=ANORM*4.D0/27.D0/DSQRT(6.D0)*
     * X*(2.D0-X/3.D0)*DEXP(-X/3.D0)
      IF(L.EQ.2)RINIT=ANORM*4.D0/81.D0/DSQRT(30.D0)*
     * X**2*DEXP(-X/3.D0)
      RETURN
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      FUNCTION RFINAL(E,L,X,Z)
C
C   HYDROGENIC RADIAL WAVE FUNCTION FOR E>0 AND E=0.
C   E: ENERGY IN ATOMIC UNITS
C   L: ANGULAR MOMENTUM QUANTUM NUMBER
C   X=Z*R, WHERE:
C        Z: ATOMIC NUMBER
C        R: RADIAL DISTANCE IN ATOMIC UNITS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 XX,ETA1,ZLMIN,UNITR
      IF(X.EQ.0.D0)X=1.D-20
      UNITR=(1.D0,0.D0)
      NL=1
      IFAIL=1
      IF(E.EQ.0.D0)GO TO 1
      PI=3.1415926535897932D0
      AK=DSQRT(2.D0*E)
      ETA1=-Z/AK*UNITR
      ZLMIN=DFLOAT(L)*UNITR
      MODE1=1
      KFN=-1
      ANORM=Z*DSQRT(2.D0/PI/AK)
      XX=AK*X/Z*UNITR
      XV=X/Z
      CALL SCOUL(-Z,E,L,XV,F,FP,G,GP,ERR)
      RFINAL=ANORM/X*F
      RETURN
1     CONTINUE
      ETA1=(0.D0,0.D0)
      ZLMIN=DFLOAT(2*L+1)*UNITR
      MODE1=22
      KFN=2
      XX=DSQRT(8.D0*X)*UNITR
      XBESS=DSQRT(8.D0*X)
      LORDER=2*L+1
      RFINAL=DSQRT(2.D0*Z/X)*BESSJ(LORDER,XBESS)
      RETURN
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       ADBP2299
CCCCCCCCCC         COULOMB AND BESSEL FUNCTIONS        CCCCCCCCCC       ADBP2300
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       ADBP2301
C  **************************************************************       ADBP2302
C                       SUBROUTINE SCOUL                                ADBP2303
C  **************************************************************       ADBP2304
      SUBROUTINE SCOUL(Z,E,L,R,F,FP,G,GP,ERR)                           ADBP2305
C                                                                       ADBP2306
C     THIS SUBROUTINE COMPUTES RADIAL SCHRODINGER-COULOMB WAVE          ADBP2307
C  FUNCTIONS FOR FREE STATES.                                           ADBP2308
C                                                                       ADBP2309
C  **** ALL QUANTITIES IN ATOMIC UNITS.                                 ADBP2310
C                                                                       ADBP2311
C  INPUT ARGUMENTS:                                                     ADBP2312
C     Z ........ FIELD STRENGTH, I.E. VALUE OF R*V(R) (ASSUMED          ADBP2313
C                CONSTANT).                                             ADBP2314
C     E ........ PARTICLE KINETIC ENERGY (POSITIVE).                    ADBP2315
C     K ........ ANGULAR MOMENTUM QUANTUM NUMBER KAPPA (.NE.0).         ADBP2316
C     R ........ RADIAL DISTANCE (POSITIVE).                            ADBP2317
C                                                                       ADBP2318
C  OUTPUT ARGUMENTS:                                                    ADBP2319
C     F, FP .... REGULAR RADIAL SCHRODINGER-COULOMB FUNCTION AND        ADBP2320
C                ITS DERIVATIVE.                                        ADBP2321
C     G, GP .... IRREGULAR RADIAL SCHRODINGER-COULOMB FUNCTION          ADBP2322
C                AND ITS DERIVATIVE.                                    ADBP2323
C     ERR ...... ACCURACY OF THE COMPUTED FUNCTIONS (RELATIVE UN-       ADBP2324
C                CERTAINTY).                                            ADBP2325
C  OUTPUT THROUGH COMMON/OCOUL/:                                        ADBP2326
C     WAVNUM ... WAVE NUMBER.                                           ADBP2327
C     ETA ...... SOMMERFELD'S PARAMETER.                                ADBP2328
C     DELTA .... COULOMB PHASE SHIFT (MODULUS 2*PI).                    ADBP2329
C                                                                       ADBP2330
C     RADIAL FUNCTIONS ARE NORMALIZED SO THAT, FOR LARGE R, THEY        ADBP2331
C  OSCILLATE WITH UNIT AMPLITUDE.                                       ADBP2332
C                                                                       ADBP2333
C     OTHER SUBPROGRAMS REQUIRED: SUBROUTINES FCOUL AND SUM2F0,         ADBP2334
C                                 AND FUNCTION CLGAM.                   ADBP2335
C                                                                       ADBP2336
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           ADBP2337
      COMMON/OFCOUL/DELTA0                                              ADBP2338
      COMMON/OCOUL/WAVNUM,ETA,DELTA                                     ADBP2339
C                                                                       ADBP2340
      IF(E.LT.0.000001D0.OR.L.LT.0) THEN                                  ADBP2341
        F=0.0D0                                                         ADBP2342
        FP=0.0D0                                                        ADBP2343
        G=1.0D35                                                        ADBP2344
        GP=-1.0D35                                                      ADBP2345
        ERR=1.0D0                                                       ADBP2346
        IF(E.LT.0.000001D0) WRITE(6,2101)                                 ADBP2347
 2101   FORMAT(1X,'*** ERROR IN SCOUL: E IS TOO SMALL.')                ADBP2348
        IF(L.LT.0) WRITE(6,2102)                                        ADBP2349
 2102   FORMAT(1X,'*** ERROR IN SCOUL: L.LT.0.')                        ADBP2350
        RETURN                                                          ADBP2351
      ENDIF                                                             ADBP2352
C                                                                       ADBP2353
C  ****  PARAMETERS.                                                    ADBP2354
C                                                                       ADBP2355
      WAVNUM=DSQRT(E+E)                                                 ADBP2356
      IF(DABS(Z).GT.0.00001D0) THEN                                     ADBP2357
        ETA=Z/WAVNUM                                                    ADBP2358
        ICAL=0                                                          ADBP2359
      ELSE                                                              ADBP2360
        ETA=0.0D0                                                       ADBP2361
        ICAL=1                                                          ADBP2362
      ENDIF                                                             ADBP2363
      RLAMB=L                                                           ADBP2364
      X=WAVNUM*R                                                        ADBP2365
      IF(ICAL.EQ.1) GO TO 1                                             ADBP2366
C                                                                       ADBP2367
C  ************  COULOMB FUNCTIONS.                                     ADBP2368
C                                                                       ADBP2369
      DELTA0=1.0D30                                                     ADBP2370
      CALL FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)                             ADBP2371
      FP=FP*WAVNUM                                                      ADBP2372
      GP=GP*WAVNUM                                                      ADBP2373
      IF(DELTA0.LT.1.0D29) THEN                                         ADBP2374
        DELTA=DELTA0                                                    ADBP2375
      ELSE                                                              ADBP2376
        DELTA=DELTAC(ETA,RLAMB)                                         ADBP2377
      ENDIF                                                             ADBP2378
      RETURN                                                            ADBP2379
C                                                                       ADBP2380
C  ************  Z=0. SPHERICAL BESSEL FUNCTIONS.                       ADBP2381
C                                                                       ADBP2382
    1 CONTINUE                                                          ADBP2383
      call besselJN(1,L,X,BESJN)
      F=X*BESJN                                                         ADBP2384
      call besselJN(2,L,X,BESJN)
      G=-X*BESJN                                                        ADBP2385

      call besselJN(1,L,X,BESJN)
      call besselJN(1,L+1,X,BESJN_L)

      FP=((L+1)*BESJN-X*BESJN_L)*WAVNUM                                 ADBP2386

      call besselJN(2,L,X,BESJN)
      call besselJN(2,L+1,X,BESJN_L)

      GP=-((L+1)*BESJN-X*BESJN_L)*WAVNUM                                ADBP2387
      DELTA=0.0D0                                                       ADBP2388
      ERR=0.0D0                                                         ADBP2389
      RETURN                                                            ADBP2390
      END                                                               ADBP2391
C  **************************************************************       ADBP2392 

C  **************************************************************       ADBP2535
C                       SUBROUTINE FCOUL                                ADBP2536
C  **************************************************************       ADBP2537
      SUBROUTINE FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)                       ADBP2538
C                                                                       ADBP2539
C     CALCULATION OF COULOMB FUNCTIONS FOR REAL ETA, RLAMB.GT.-1        ADBP2540
C  AND X LARGER THAN, OR OF THE ORDER OF XTP0 (THE TURNING POINT        ADBP2541
C  FOR RLAMB=0). STEED'S CONTINUED FRACTION METHOD IS COMBINED          ADBP2542
C  WITH RECURSION RELATIONS AND AN ASYMPTOTIC EXPANSION. THE            ADBP2543
C  OUTPUT VALUE ERR=1.0D0 INDICATES THAT THE ADOPTED EVALUATION         ADBP2544
C  ALGORITHM IS NOT APPLICABLE (X IS TOO SMALL).                        ADBP2545
C                                                                       ADBP2546
C  INPUT ARGUMENTS:                                                     ADBP2547
C     ETA ...... SOMMERFELD'S PARAMETER.                                ADBP2548
C     RLAMB .... ANGULAR MOMENTUM.                                      ADBP2549
C     X ........ VARIABLE (=WAVE NUMBER TIMES RADIAL DISTANCE).         ADBP2550
C                                                                       ADBP2551
C  OUTPUT ARGUMENTS:                                                    ADBP2552
C     F, FP .... REGULAR FUNCTION AND ITS DERIVATIVE.                   ADBP2553
C     G, GP .... IRREGULAR FUNCTION AND ITS DERIVATIVE.                 ADBP2554
C     ERR ...... RELATIVE NUMERICAL UNCERTAINTY. A VALUE OF THE         ADBP2555
C                ORDER OF 10**(-N) MEANS THAT THE CALCULATED            ADBP2556
C                FUNCTIONS ARE ACCURATE TO N DECIMAL FIGURES.           ADBP2557
C                THE MAXIMUM ACCURACY ATTAINABLE WITH DOUBLE            ADBP2558
C                PRECISION ARITHMETIC IS ABOUT 1.0D-15.                 ADBP2559
C                                                                       ADBP2560
C     OTHER SUBPROGRAMS REQUIRED: SUBROUTINE SUM2F0 AND                 ADBP2561
C                                 FUNCTIONS DELTAC AND CLGAM.           ADBP2562
C                                                                       ADBP2563
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           ADBP2564
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI,TPI=PI+PI,        ADBP2565
     1  EPS=1.0D-16,TOP=1.0D5,NTERM=1000)                               ADBP2566
      COMMON/OFCOUL/DELTA                                               ADBP2567
C                                                                       ADBP2568
      IF(RLAMB.LT.-0.999D0) THEN                                        ADBP2569
        WRITE(6,'(1X,''*** ERROR IN RCOUL: RLAMB.LT.-0.999'')')         ADBP2570
        STOP                                                            ADBP2571
      ENDIF                                                             ADBP2572
      IF(X.LT.EPS) GO TO 10                                             ADBP2573
C                                                                       ADBP2574
C  ****  NUMERICAL CONSTANTS.                                           ADBP2575
C                                                                       ADBP2576
      CI=DCMPLX(0.0D0,1.0D0)                                            ADBP2577
      CI2=2.0D0*CI                                                      ADBP2578
      CIETA=CI*ETA                                                      ADBP2579
      X2=X*X                                                            ADBP2580
      ETA2=ETA*ETA                                                      ADBP2581
C                                                                       ADBP2582
C  ****  TURNING POINT (XTP). (44)                                      ADBP2583
C                                                                       ADBP2584
      IF(RLAMB.GE.0.0D0) THEN                                           ADBP2585
        XTP=ETA+DSQRT(ETA2+RLAMB*(RLAMB+1.0D0))                         ADBP2586
      ELSE                                                              ADBP2587
        XTP=EPS                                                         ADBP2588
      ENDIF                                                             ADBP2589
      ERRS=10.0D0                                                       ADBP2590
      IF(X.LT.XTP) GO TO 1                                              ADBP2591
C                                                                       ADBP2592
C  ************  ASYMPTOTIC EXPANSION. (71-75)                          ADBP2593
C                                                                       ADBP2594
C  ****  COULOMB PHASE-SHIFT.                                           ADBP2595
      DELTA=DELTAC(ETA,RLAMB)                                           ADBP2596
C                                                                       ADBP2597
      CPA=CIETA-RLAMB                                                   ADBP2598
      CPB=CIETA+RLAMB+1.0D0                                             ADBP2599
      CPZ=CI2*X                                                         ADBP2600
      CALL SUM2F0(CPA,CPB,CPZ,C2F0,ERR1)                                ADBP2601
      CQA=CPA+1.0D0                                                     ADBP2602
      CQB=CPB+1.0D0                                                     ADBP2603
      CALL SUM2F0(CQA,CQB,CPZ,C2F0P,ERR2)                               ADBP2604
      C2F0P=CI*C2F0P*CPA*CPB/(2.0D0*X2)                                 ADBP2605
C  ****  FUNCTIONS.                                                     ADBP2606
      THETA=X-ETA*DLOG(2.0D0*X)-RLAMB*PIH+DELTA                         ADBP2607
      IF(THETA.GT.1.0D4) THETA=DMOD(THETA,TPI)                          ADBP2608
      CEITH=CDEXP(CI*THETA)                                             ADBP2609
      CGIF=C2F0*CEITH                                                   ADBP2610
      G=CGIF                                                            ADBP2611
      F=-CI*CGIF                                                        ADBP2612
C  ****  DERIVATIVES.                                                   ADBP2613
      CGIFP=(C2F0P+CI*(1.0D0-ETA/X)*C2F0)*CEITH                         ADBP2614
      GP=CGIFP                                                          ADBP2615
      FP=-CI*CGIFP                                                      ADBP2616
C  ****  GLOBAL UNCERTAINTY. THE WRONSKIAN MAY DIFFER FROM 1 DUE        ADBP2617
C        TO TRUNCATION AND ROUNDOFF ERRORS.                             ADBP2618
      ERR=DMAX1(ERR1,ERR2,DABS(G*FP-F*GP-1.0D0))                        ADBP2619
      IF(ERR.LE.EPS) RETURN                                             ADBP2620
      ERRS=ERR                                                          ADBP2621
C                                                                       ADBP2622
C  ************  STEED'S CONTINUED FRACTION METHOD.                     ADBP2623
C                                                                       ADBP2624
    1 CONTINUE                                                          ADBP2625
      CIETA2=CIETA+CIETA                                                ADBP2626
      ETAX=ETA*X                                                        ADBP2627
C                                                                       ADBP2628
C  ****  CONTINUED FRACTION FOR F. (60-70)                              ADBP2629
C                                                                       ADBP2630
      INULL=0                                                           ADBP2631
      RLAMBN=RLAMB+1.0D0                                                ADBP2632
      A1=-(RLAMBN+1.0D0)*(RLAMBN**2+ETA2)*X/RLAMBN                      ADBP2633
      B0=(RLAMBN/X)+(ETA/RLAMBN)                                        ADBP2634
      B1=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)              ADBP2635
      FA3=B0                                                            ADBP2636
      FA2=B0*B1+A1                                                      ADBP2637
      FB3=1.0D0                                                         ADBP2638
      FB2=B1                                                            ADBP2639
      RF=FA3                                                            ADBP2640
C                                                                       ADBP2641
      DO 2 N=2,NTERM                                                    ADBP2642
      RFO=RF                                                            ADBP2643
      DAF=DABS(RF)                                                      ADBP2644
      RLAMBN=RLAMB+N                                                    ADBP2645
      AN=-(RLAMBN**2-1.0D0)*(RLAMBN**2+ETA2)*X2                         ADBP2646
      BN=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)              ADBP2647
      FA1=FA2*BN+FA3*AN                                                 ADBP2648
      FB1=FB2*BN+FB3*AN                                                 ADBP2649
      TST=DABS(FB1)                                                     ADBP2650
C                                                                       ADBP2651
      IF(TST.LT.1.0D-25) THEN                                           ADBP2652
        IF(INULL.GT.0) STOP                                             ADBP2653
        INULL=1                                                         ADBP2654
        FA3=FA2                                                         ADBP2655
        FA2=FA1                                                         ADBP2656
        FB3=FB2                                                         ADBP2657
        FB2=FB1                                                         ADBP2658
        RF=RFO                                                          ADBP2659
      ELSE                                                              ADBP2660
        FA3=FA2/TST                                                     ADBP2661
        FA2=FA1/TST                                                     ADBP2662
        FB3=FB2/TST                                                     ADBP2663
        FB2=FB1/TST                                                     ADBP2664
        RF=FA2/FB2                                                      ADBP2665
        IF(DABS(RF-RFO).LT.EPS*DAF) GO TO 3                             ADBP2666
      ENDIF                                                             ADBP2667
    2 CONTINUE                                                          ADBP2668
    3 CONTINUE                                                          ADBP2669
      IF(DAF.GT.1.0D-25) THEN                                           ADBP2670
        ERRF=DABS(RF-RFO)/DAF                                           ADBP2671
      ELSE                                                              ADBP2672
        ERRF=EPS                                                        ADBP2673
      ENDIF                                                             ADBP2674
      IF(ERRF.GT.ERRS) THEN                                             ADBP2675
        ERR=ERRS                                                        ADBP2676
        RETURN                                                          ADBP2677
      ENDIF                                                             ADBP2678
C                                                                       ADBP2679
C  ****  DOWNWARD RECURSION FOR F AND FP. ONLY IF RLAMB.GT.1 AND        ADBP2680
C        X.LT.XTP. (48,49)                                              ADBP2681
C                                                                       ADBP2682
      RLAMB0=RLAMB                                                      ADBP2683
      IF(X.GE.XTP.OR.RLAMB0.LT.1.0D0) THEN                              ADBP2684
        ISHIFT=0                                                        ADBP2685
      ELSE                                                              ADBP2686
        FT=1.0D0                                                        ADBP2687
        FTP=RF                                                          ADBP2688
        IS0=RLAMB0+1.0D-6                                               ADBP2689
        TST=X*(X-2.0D0*ETA)                                             ADBP2690
        DO 4 I=1,IS0                                                    ADBP2691
        ETARL0=ETA/RLAMB0                                               ADBP2692
        RL=DSQRT(1.0D0+ETARL0**2)                                       ADBP2693
        SL=(RLAMB0/X)+ETARL0                                            ADBP2694
        RLAMB0=RLAMB0-1.0D0                                             ADBP2695
        FTO=FT                                                          ADBP2696
        FT=(SL*FT+FTP)/RL                                               ADBP2697
        FTP=SL*FT-RL*FTO                                                ADBP2698
        IF(FT.GT.1.0D10) THEN                                           ADBP2699
          FTP=FTP/FT                                                    ADBP2700
          FT=1.0D0                                                      ADBP2701
        ENDIF                                                           ADBP2702
        RL1T=RLAMB0*(RLAMB0+1.0D0)                                      ADBP2703
        IF(TST.GT.RL1T) THEN                                            ADBP2704
          ISHIFT=I                                                      ADBP2705
          GO TO 5                                                       ADBP2706
        ENDIF                                                           ADBP2707
    4   CONTINUE                                                        ADBP2708
        ISHIFT=IS0                                                      ADBP2709
    5   CONTINUE                                                        ADBP2710
        XTPC=ETA+DSQRT(ETA2+RL1T)                                       ADBP2711
        RFM=FTP/FT                                                      ADBP2712
      ENDIF                                                             ADBP2713
C                                                                       ADBP2714
C  ****  CONTINUED FRACTION FOR P+CI*Q WITH RLAMB0. (76-79)             ADBP2715
C                                                                       ADBP2716
      INULL=0                                                           ADBP2717
      CAN=CIETA-ETA2-RLAMB0*(RLAMB0+1.0D0)                              ADBP2718
      CB0=X-ETA                                                         ADBP2719
      CBN=2.0D0*(X-ETA+CI)                                              ADBP2720
      CFA3=CB0                                                          ADBP2721
      CFA2=CB0*CBN+CAN                                                  ADBP2722
      CFB3=1.0D0                                                        ADBP2723
      CFB2=CBN                                                          ADBP2724
      CPIQ=CFA3                                                         ADBP2725
C                                                                       ADBP2726
      DO 6 N=2,NTERM                                                    ADBP2727
      CPIQO=CPIQ                                                        ADBP2728
      DAPIQ=CDABS(CPIQ)                                                 ADBP2729
      CAN=CAN+CIETA2+(N+N-2)                                            ADBP2730
      CBN=CBN+CI2                                                       ADBP2731
      CFA1=CFA2*CBN+CFA3*CAN                                            ADBP2732
      CFB1=CFB2*CBN+CFB3*CAN                                            ADBP2733
      TST=CDABS(CFB1)                                                   ADBP2734
C                                                                       ADBP2735
      IF(TST.LT.1.0D-25) THEN                                           ADBP2736
        IF(INULL.GT.0) STOP                                             ADBP2737
        INULL=1                                                         ADBP2738
        CFA3=CFA2                                                       ADBP2739
        CFA2=CFA1                                                       ADBP2740
        CFB3=CFB2                                                       ADBP2741
        CFB2=CFB1                                                       ADBP2742
        CPIQ=CPIQO                                                      ADBP2743
      ELSE                                                              ADBP2744
        CFA3=CFA2/TST                                                   ADBP2745
        CFA2=CFA1/TST                                                   ADBP2746
        CFB3=CFB2/TST                                                   ADBP2747
        CFB2=CFB1/TST                                                   ADBP2748
        CPIQ=CFA2/CFB2                                                  ADBP2749
        IF(CDABS(CPIQ-CPIQO).LT.EPS*DAPIQ) GO TO 7                      ADBP2750
      ENDIF                                                             ADBP2751
    6 CONTINUE                                                          ADBP2752
    7 CONTINUE                                                          ADBP2753
      IF(DAPIQ.GT.1.0D-25) THEN                                         ADBP2754
        ERRPIQ=CDABS(CPIQ-CPIQO)/DAPIQ                                  ADBP2755
      ELSE                                                              ADBP2756
        ERRPIQ=EPS                                                      ADBP2757
      ENDIF                                                             ADBP2758
      IF(ERRPIQ.GT.ERRS) THEN                                           ADBP2759
        ERR=ERRS                                                        ADBP2760
        RETURN                                                          ADBP2761
      ENDIF                                                             ADBP2762
      CPIQ=CI*CPIQ/X                                                    ADBP2763
C                                                                       ADBP2764
      RP=CPIQ                                                           ADBP2765
      RQ=-CI*CPIQ                                                       ADBP2766
      IF(RQ.LE.1.0D-25) GO TO 10                                        ADBP2767
      ERR=DMAX1(ERRF,ERRPIQ)                                            ADBP2768
C                                                                       ADBP2769
C  ****  INVERTING STEED'S TRANSFORMATION. (57,58)                      ADBP2770
C                                                                       ADBP2771
      IF(ISHIFT.LT.1) THEN                                              ADBP2772
        RFP=RF-RP                                                       ADBP2773
        F=DSQRT(RQ/(RFP**2+RQ**2))                                      ADBP2774
        IF(FB2.LT.0.0D0) F=-F                                           ADBP2775
        FP=RF*F                                                         ADBP2776
        G=RFP*F/RQ                                                      ADBP2777
        GP=(RP*RFP-RQ**2)*F/RQ                                          ADBP2778
        IF(X.LT.XTP.AND.G.GT.TOP*F) GO TO 10                            ADBP2779
      ELSE                                                              ADBP2780
        RFP=RFM-RP                                                      ADBP2781
        FM=DSQRT(RQ/(RFP**2+RQ**2))                                     ADBP2782
        G=RFP*FM/RQ                                                     ADBP2783
        GP=(RP*RFP-RQ**2)*FM/RQ                                         ADBP2784
        IF(X.LT.XTPC.AND.G.GT.TOP*FM) GO TO 10                          ADBP2785
C  ****  UPWARD RECURSION FOR G AND GP (IF ISHIFT.GT.0). (50,51)        ADBP2786
        DO 8 I=1,ISHIFT                                                 ADBP2787
        RLAMB0=RLAMB0+1.0D0                                             ADBP2788
        ETARL0=ETA/RLAMB0                                               ADBP2789
        RL=DSQRT(1.0D0+ETARL0**2)                                       ADBP2790
        SL=(RLAMB0/X)+ETARL0                                            ADBP2791
        GO=G                                                            ADBP2792
        G=(SL*GO-GP)/RL                                                 ADBP2793
        GP=RL*GO-SL*G                                                   ADBP2794
        IF(G.GT.1.0D35) GO TO 10                                        ADBP2795
    8   CONTINUE                                                        ADBP2796
    9   W=RF*G-GP                                                       ADBP2797
        F=1.0D0/W                                                       ADBP2798
        FP=RF/W                                                         ADBP2799
      ENDIF                                                             ADBP2800
C  ****  THE WRONSKIAN MAY DIFFER FROM 1 DUE TO ROUNDOFF ERRORS.        ADBP2801
      ERR=DMAX1(ERR,DABS(FP*G-F*GP-1.0D0))                              ADBP2802
      RETURN                                                            ADBP2803
C                                                                       ADBP2804
   10 F=0.0D0                                                           ADBP2805
      FP=0.0D0                                                          ADBP2806
      G=1.0D35                                                          ADBP2807
      GP=-1.0D35                                                        ADBP2808
      ERR=1.0D0                                                         ADBP2809
      RETURN                                                            ADBP2810
      END                                                               ADBP2811
C  **************************************************************       ADBP2812
C                       SUBROUTINE SUM2F0                               ADBP2813
C  **************************************************************       ADBP2814
      SUBROUTINE SUM2F0(CA,CB,CZ,CF,ERR)                                ADBP2815
C                                                                       ADBP2816
C     SUMMATION OF THE 2F0(CA,CB;CS) HYPERGEOMETRIC ASYMPTOTIC          ADBP2817
C  SERIES. THE POSITIVE AND NEGATIVE CONTRIBUTIONS TO THE REAL          ADBP2818
C  AND IMAGINARY PARTS ARE ADDED SEPARATELY TO OBTAIN AN ESTIMATE       ADBP2819
C  OF ROUNDING ERRORS.                                                  ADBP2820
C                                                                       ADBP2821
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           ADBP2822
      PARAMETER (EPS=1.0D-16,ACCUR=0.5D-15,NTERM=75)                    ADBP2823
      RRP=1.0D0                                                         ADBP2824
      RRN=0.0D0                                                         ADBP2825
      RIP=0.0D0                                                         ADBP2826
      RIN=0.0D0                                                         ADBP2827
      CDF=1.0D0                                                         ADBP2828
      ERR2=0.0D0                                                        ADBP2829
      ERR3=1.0D0                                                        ADBP2830
      DO 1 I=1,NTERM                                                    ADBP2831
      J=I-1                                                             ADBP2832
      CDF=CDF*(CA+J)*(CB+J)/(I*CZ)                                      ADBP2833
      ERR1=ERR2                                                         ADBP2834
      ERR2=ERR3                                                         ADBP2835
      ERR3=CDABS(CDF)                                                   ADBP2836
      IF(ERR1.GT.ERR2.AND.ERR2.LT.ERR3) GO TO 2                         ADBP2837
      AR=CDF                                                            ADBP2838
      IF(AR.GT.0.0D0) THEN                                              ADBP2839
        RRP=RRP+AR                                                      ADBP2840
      ELSE                                                              ADBP2841
        RRN=RRN+AR                                                      ADBP2842
      ENDIF                                                             ADBP2843
      AI=DCMPLX(0.0D0,-1.0D0)*CDF                                       ADBP2844
      IF(AI.GT.0.0D0) THEN                                              ADBP2845
        RIP=RIP+AI                                                      ADBP2846
      ELSE                                                              ADBP2847
        RIN=RIN+AI                                                      ADBP2848
      ENDIF                                                             ADBP2849
      CF=DCMPLX(RRP+RRN,RIP+RIN)                                        ADBP2850
      AF=CDABS(CF)                                                      ADBP2851
      IF(AF.GT.1.0D25) THEN                                             ADBP2852
        CF=0.0D0                                                        ADBP2853
        ERR=1.0D0                                                       ADBP2854
        RETURN                                                          ADBP2855
      ENDIF                                                             ADBP2856
      IF(ERR3.LT.1.0D-25*AF.OR.ERR3.LT.EPS) THEN                        ADBP2857
         ERR=EPS                                                        ADBP2858
         RETURN                                                         ADBP2859
      ENDIF                                                             ADBP2860
    1 CONTINUE                                                          ADBP2861
C  ****  ROUNDOFF ERROR.                                                ADBP2862
    2 CONTINUE                                                          ADBP2863
      TR=DABS(RRP+RRN)                                                  ADBP2864
      IF(TR.GT.1.0D-25) THEN                                            ADBP2865
        ERRR=(RRP-RRN)*ACCUR/TR                                         ADBP2866
      ELSE                                                              ADBP2867
        ERRR=1.0D0                                                      ADBP2868
      ENDIF                                                             ADBP2869
      TI=DABS(RIP+RIN)                                                  ADBP2870
      IF(TI.GT.1.0D-25) THEN                                            ADBP2871
        ERRI=(RIP-RIN)*ACCUR/TI                                         ADBP2872
      ELSE                                                              ADBP2873
        ERRI=1.0D0                                                      ADBP2874
      ENDIF                                                             ADBP2875
C  ****  ... AND TRUNCATION ERROR.                                      ADBP2876
      IF(AR.GT.1.0D-25) THEN                                            ADBP2877
      ERR=DMAX1(ERRR,ERRI)+ERR2/AF                                      ADBP2878
      ELSE                                                              ADBP2879
      ERR=DMAX1(ERRR,ERRI)                                              ADBP2880
      ENDIF                                                             ADBP2881
      RETURN                                                            ADBP2882
      END                                                               ADBP2883
C  **************************************************************       ADBP2884
C                         FUNCTION DELTAC                               ADBP2885
C  **************************************************************       ADBP2886
      FUNCTION DELTAC(ETA,RLAMB)                                        ADBP2887
C                                                                       ADBP2888
C     CALCULATION OF COULOMB PHASE SHIFT (MODULUS 2*PI). (47)           ADBP2889
C                                                                       ADBP2890
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           ADBP2891
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)                     ADBP2892
      CI=DCMPLX(0.0D0,1.0D0)                                            ADBP2893
C  ****  COULOMB PHASE-SHIFT.                                           ADBP2894
      DELTAC=-CI*CLGAM(RLAMB+1.0D0+CI*ETA)                              ADBP2895
      IF(DELTAC.GE.0.0D0) THEN                                          ADBP2896
        DELTAC=DMOD(DELTAC,TPI)                                         ADBP2897
      ELSE                                                              ADBP2898
        DELTAC=-DMOD(-DELTAC,TPI)                                       ADBP2899
      ENDIF                                                             ADBP2900
      RETURN                                                            ADBP2901
      END                                                               ADBP2902
C  **************************************************************       ADBP2903
C                       FUNCTION CLGAM                                  ADBP2904
C  **************************************************************       ADBP2905
      FUNCTION CLGAM(CZ)                                                ADBP2906
C                                                                       ADBP2907
C     THIS FUNCTION GIVES LOG(GAMMA(CZ)) FOR COMPLEX ARGUMENTS.         ADBP2908
C                                                                       ADBP2909
C   REF.: M. ABRAMOWITZ AND I.A. STEGUN, 'HANDBOOK OF MATHEMATI-        ADBP2910
C         CAL FUNCTIONS'. DOVER, NEW YORK (1974). PP 255-257.           ADBP2911
C                                                                       ADBP2912
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           ADBP2913
      PARAMETER (PI=3.1415926535897932D0)                               ADBP2914
      CZA=CZ                                                            ADBP2915
      ICONJ=0                                                           ADBP2916
      AR=CZA                                                            ADBP2917
      CLGAM=36.84136149D0                                               ADBP2918
      IF(CDABS(CZA).LT.1.0D-16) RETURN                                  ADBP2919
C                                                                       ADBP2920
      AI=CZA*DCMPLX(0.0D0,-1.0D0)                                       ADBP2921
      IF(AI.GT.0.0D0) THEN                                              ADBP2922
        ICONJ=0                                                         ADBP2923
      ELSE                                                              ADBP2924
        ICONJ=1                                                         ADBP2925
        CZA=DCONJG(CZA)                                                 ADBP2926
      ENDIF                                                             ADBP2927
C                                                                       ADBP2928
      CZFAC=1.0D0                                                       ADBP2929
      CZFL=0.0D0                                                        ADBP2930
    1 CZFAC=CZFAC/CZA                                                   ADBP2931
      IF(CDABS(CZFAC).GT.1.0D8) THEN                                    ADBP2932
        CZFL=CZFL+CDLOG(CZFAC)                                          ADBP2933
        CZFAC=1.0D0                                                     ADBP2934
      ENDIF                                                             ADBP2935
      CZA=CZA+1.0D0                                                     ADBP2936
      AR=CZA                                                            ADBP2937
      IF(CDABS(CZA).LT.1.0D-16) RETURN                                  ADBP2938
      IF(CDABS(CZA).GT.15.0D0.AND.AR.GT.0.0D0) GO TO 2                  ADBP2939
      GO TO 1                                                           ADBP2940
C  ****  STIRLING'S EXPANSION OF CDLOG(GAMMA(CZA)).                     ADBP2941
    2 CZI2=1.0D0/(CZA*CZA)                                              ADBP2942
      CZS=(43867.0D0/244188.0D0)*CZI2                                   ADBP2943
      CZS=(CZS-3617.0D0/122400.0D0)*CZI2                                ADBP2944
      CZS=(CZS+1.0D0/156.0D0)*CZI2                                      ADBP2945
      CZS=(CZS-691.0D0/360360.0D0)*CZI2                                 ADBP2946
      CZS=(CZS+1.0D0/1188.0D0)*CZI2                                     ADBP2947
      CZS=(CZS-1.0D0/1680.0D0)*CZI2                                     ADBP2948
      CZS=(CZS+1.0D0/1260.0D0)*CZI2                                     ADBP2949
      CZS=(CZS-1.0D0/360.0D0)*CZI2                                      ADBP2950
      CZS=(CZS+1.0D0/12.0D0)/CZA                                        ADBP2951
      CLGAM=(CZA-0.5D0)*CDLOG(CZA)-CZA+9.1893853320467274D-1+CZS        ADBP2952
     1     +CZFL+CDLOG(CZFAC)                                           ADBP2953
      IF(ICONJ.EQ.1) CLGAM=DCONJG(CLGAM)                                ADBP2954
      RETURN                                                            ADBP2955
      END                                                               ADBP2956
C  **************************************************************       ADBP2957
C  **************************************************************       ADBP2957
C                         FUNCION BESJN                                 ADBP2958
C  **************************************************************       ADBP2959
       subroutine besselJN(JY,N,X,BESJN)                                ADBP2960
C                                                                       ADBP2961
C      THIS FUNCTION COMPUTES THE SPHERICAL BESSEL FUNCTIONS OF         ADBP2962
C   THE FIRST KIND AND SPHERICAL BESSEL FUNCTIONS OF THE SECOND         ADBP2963
C   KIND (ALSO KNOWN AS SPHERICAL NEUMANN FUNCTIONS) FOR REAL           ADBP2964
C   POSITIVE ARGUMENTS.                                                 ADBP2965
C                                                                       ADBP2966
C      INPUT:                                                           ADBP2967
C         JY ...... KIND: 1(BESSEL) OR 2(NEUMANN).                      ADBP2968
C         N ....... ORDER (INTEGER).                                    ADBP2969
C         X ....... ARGUMENT (REAL AND POSITIVE).                       ADBP2970
C                                                                       ADBP2971
C   REF.: M. ABRAMOWITZ AND I.A. STEGUN, 'HANDBOOK OF MATHEMATI-        ADBP2972
C         CAL FUNCTIONS'. DOVER, NEW YORK (1974). PP 435-478.           ADBP2973
C                                                                       ADBP2974
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               ADBP2975
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)                     ADBP2976
      IF(X.LT.0) THEN                                                   ADBP2977
        WRITE(6,1000)                                                   ADBP2978
 1000   FORMAT(1X,'*** NEGATIVE ARGUMENT IN FUNCTION BESJN.')           ADBP2979
        STOP                                                            ADBP2980
      ENDIF                                                             ADBP2981
C  ****  ORDER AND PHASE CORRECTION FOR NEUMANN FUNCTIONS.              ADBP2982
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.15.                            ADBP2983
      IF(JY.EQ.2) THEN                                                  ADBP2984
        NL=-N-1                                                         ADBP2985
        IPH=2*MOD(IABS(N),2)-1                                          ADBP2986
      ELSE                                                              ADBP2987
        NL=N                                                            ADBP2988
        IPH=1                                                           ADBP2989
      ENDIF                                                             ADBP2990
C  ****  SELECTION OF CALCULATION MODE.                                 ADBP2991
      IF(NL.LT.0) GO TO 10                                              ADBP2992
      IF(X.GT.1.0D0*NL) GO TO 7                                         ADBP2993
      XI=X*X                                                            ADBP2994
      IF(XI.GT.NL+NL+3.0D0) GO TO 4                                     ADBP2995
C  ****  POWER SERIES FOR SMALL ARGUMENTS AND POSITIVE ORDERS.          ADBP2996
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.2.                             ADBP2997
      F1=1.0D0                                                          ADBP2998
      IP=1                                                              ADBP2999
      IF(NL.NE.0) THEN                                                  ADBP3000
        DO 1 I=1,NL                                                     ADBP3001
        IP=IP+2                                                         ADBP3002
    1   F1=F1*X/IP                                                      ADBP3003
      ENDIF                                                             ADBP3004
      XI=0.5D0*XI                                                       ADBP3005
      BESJN=1.0D0                                                       ADBP3006
      PS=1.0D0                                                          ADBP3007
      DO 2 I=1,500                                                      ADBP3008
      IP=IP+2                                                           ADBP3009
      PS=-PS*XI/(I*IP)                                                  ADBP3010
      BESJN=BESJN+PS                                                    ADBP3011
      IF(DABS(PS).LT.1.0D-18*DABS(BESJN)) GO TO 3                       ADBP3012
    2 CONTINUE                                                          ADBP3013
    3 BESJN=IPH*F1*BESJN                                                ADBP3014
      RETURN                                                            ADBP3015
C  ****  MILLER'S METHOD FOR POSITIVE ORDERS AND INTERMEDIATE           ADBP3016
C        ARGUMENTS. ABRAMOWITZ AND STEGUN, EQ. 10.1.19.                 ADBP3017
    4 XI=1.0D0/X                                                        ADBP3018
      F2=0.0D0                                                          ADBP3019
      F3=1.0D-35                                                        ADBP3020
      IP=2*(NL+31)+3                                                    ADBP3021
      DO 5 I=1,31                                                       ADBP3022
      F1=F2                                                             ADBP3023
      F2=F3                                                             ADBP3024
      IP=IP-2                                                           ADBP3025
      F3=IP*XI*F2-F1                                                    ADBP3026
      IF(DABS(F3).GT.1.0D30) THEN                                       ADBP3027
        F2=F2/F3                                                        ADBP3028
        F3=1.0D0                                                        ADBP3029
      ENDIF                                                             ADBP3030
    5 CONTINUE                                                          ADBP3031
      BESJN=1.0D0                                                       ADBP3032
      F2=F2/F3                                                          ADBP3033
      F3=1.0D0                                                          ADBP3034
      DO 6 I=1,NL                                                       ADBP3035
      F1=F2                                                             ADBP3036
      F2=F3                                                             ADBP3037
      IP=IP-2                                                           ADBP3038
      F3=IP*XI*F2-F1                                                    ADBP3039
      IF(DABS(F3).GT.1.0D30) THEN                                       ADBP3040
        BESJN=BESJN/F3                                                  ADBP3041
        F2=F2/F3                                                        ADBP3042
        F3=1.0D0                                                        ADBP3043
      ENDIF                                                             ADBP3044
    6 CONTINUE                                                          ADBP3045
      BESJN=IPH*XI*DSIN(X)*BESJN/F3                                     ADBP3046
      RETURN                                                            ADBP3047
C  ****  RECURRENCE RELATION FOR ARGUMENTS GREATER THAN ORDER.          ADBP3048
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.19.                            ADBP3049
    7 XI=1.0D0/X                                                        ADBP3050
      F3=XI*DSIN(X)                                                     ADBP3051
      IF(NL.EQ.0) GO TO 9                                               ADBP3052
      F2=F3                                                             ADBP3053
      F3=XI*(F2-DCOS(X))                                                ADBP3054
      IF(NL.EQ.1) GO TO 9                                               ADBP3055
      IP=1                                                              ADBP3056
      DO 8 I=2,NL                                                       ADBP3057
      F1=F2                                                             ADBP3058
      F2=F3                                                             ADBP3059
      IP=IP+2                                                           ADBP3060
    8 F3=IP*XI*F2-F1                                                    ADBP3061
    9 BESJN=IPH*F3                                                      ADBP3062
      RETURN                                                            ADBP3063
C  ****  RECURRENCE RELATION FOR NEGATIVE ORDERS.                       ADBP3064
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.19.                            ADBP3065
   10 NL=IABS(NL)                                                       ADBP3066
      IF(X.LT.7.36D-1*(NL+1)*1.0D-35**(1.0D0/(NL+1))) THEN              ADBP3067
        BESJN=-1.0D35                                                   ADBP3068
        RETURN                                                          ADBP3069
      ENDIF                                                             ADBP3070
      XI=1.0D0/X                                                        ADBP3071
      F3=XI*DSIN(X)                                                     ADBP3072
      F2=XI*(F3-DCOS(X))                                                ADBP3073
      IP=3                                                              ADBP3074
      DO 11 I=1,NL                                                      ADBP3075
      F1=F2                                                             ADBP3076
      F2=F3                                                             ADBP3077
      IP=IP-2                                                           ADBP3078
      F3=IP*XI*F2-F1                                                    ADBP3079
      IF(DABS(F3).GT.1.0D35) THEN                                       ADBP3080
        BESJN=-1.0D35                                                   ADBP3081
        RETURN                                                          ADBP3082
      ENDIF                                                             ADBP3083
   11 CONTINUE                                                          ADBP3084
      BESJN=IPH*F3                                                      ADBP3085
      RETURN                                                            ADBP3086
      END                                                               ADBP3087 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      FUNCTION BESSJ(N,X)
C      
C     THE BESSEL FUNCTION OF THE FIRST KIND OF ORDER N AND REAL       
C     ARGUMENT X
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PI=2.D0*DASIN(1.D0)
      IF(X.EQ.0.D0)X=1.D-20
      ETA=0.D0
      RLAMB=DFLOAT(N)-0.5D0
      ERR=1.0D-15
      CALL FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)
      BESSJ=0.D0
      IF(F.EQ.0.D0)RETURN
      BESSJ=DSQRT(2.D0/PI/X)*F
      RETURN
      END	
