

auto eloss_pip = [&](double pion_p, double pip_theta, double pion_det, bool outbending){
    // momentum loss correction for low momentum pions:
    // input: p = pion momentum in GeV, pip_theta = pion theta in degree, 
    //        pion_det = pion detector (2 = FD, 3 = CD),  outbending = torus polarity
    // output: dp_pion = generated momentum - reconstructed momentum = momentum loss (+) / gain (-)

    double dp_pion = 0.0;

    if(outbending == false){ // INBENDING
        if(pion_det == 2){   // Forward Detector
            if(pip_theta < 27){                                       dp_pion =  0.00342646 + (-0.00282934) *pion_p + (0.00205983)   *pow(pion_p,2) + (-0.00043158)  *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta < 27 && pion_p >= 2.5){                      dp_pion =  0.00342646 + (-0.00282934) *2.5    + (0.00205983)   *pow(2.5,2)    + (-0.00043158)  *pow(2.5,3)    + (0) *pow(2.5,4);}
            if(pip_theta > 27 && pip_theta < 28){                     dp_pion =  0.00328565 + (-0.00376042) *pion_p + (0.00433886)   *pow(pion_p,2) + (-0.00141614)  *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 27 && pip_theta < 28 && pion_p >= 1.83){   dp_pion =  0.00328565 + (-0.00376042) *1.83   + (0.00433886)   *pow(1.83,2)   + (-0.00141614)  *pow(1.83,3)   + (0) *pow(1.83,4);}
            if(pip_theta > 28 && pip_theta < 29){                     dp_pion =  0.00328579 + (-0.00281121) *pion_p + (0.00342749)   *pow(pion_p,2) + (-0.000932614) *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 28 && pip_theta < 29 && pion_p >= 2){      dp_pion =  0.00328579 + (-0.00281121) *2      + (0.00342749)   *pow(2,2)      + (-0.000932614) *pow(2,3)      + (0) *pow(2,4);}
            if(pip_theta > 29 && pip_theta < 30){                     dp_pion =  0.00167358 + (0.00441871)  *pion_p + (-0.000834667) *pow(pion_p,2) + (-0.000137968) *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 29 && pip_theta < 30 && pion_p >= 1.9){    dp_pion =  0.00167358 + (0.00441871)  *1.9    + (-0.000834667) *pow(1.9,2)    + (-0.000137968) *pow(1.9,3)    + (0) *pow(1.9,4);}
            if(pip_theta > 30 && pip_theta < 31){                     dp_pion =  0.00274159 + (0.00635686)  *pion_p + (-0.00380977)  *pow(pion_p,2) + (0.00071627)   *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 30 && pip_theta < 31 && pion_p >= 1.9){    dp_pion =  0.00274159 + (0.00635686)  *1.9    + (-0.00380977)  *pow(1.9,2)    + (0.00071627)   *pow(1.9,3)    + (0) *pow(1.9,4);}
            if(pip_theta > 31 && pip_theta < 32){                     dp_pion =  0.00450241 + (0.00248969)  *pion_p + (-0.00336795)  *pow(pion_p,2) + (0.00111193)   *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 31 && pip_theta < 32 && pion_p >= 1.8){    dp_pion =  0.00450241 + (0.00248969)  *1.8    + (-0.00336795)  *pow(1.8,2)    + (0.00111193)   *pow(1.8,3)    + (0) *pow(1.8,4);}
            if(pip_theta > 32 && pip_theta < 33){                     dp_pion =  0.00505593 + (-0.00246203) *pion_p + (0.00172984)   *pow(pion_p,2) + (-0.000406701) *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 32 && pip_theta < 33 && pion_p >= 1.8){    dp_pion =  0.00505593 + (-0.00246203) *1.8    + (0.00172984)   *pow(1.8,2)    + (-0.000406701) *pow(1.8,3)    + (0) *pow(1.8,4);}
            if(pip_theta > 33 && pip_theta < 34){                     dp_pion =  0.00273402 + (0.00440449)  *pion_p + (-0.00373488)  *pow(pion_p,2) + (0.000996612)  *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 33 && pip_theta < 34 && pion_p >= 1.8){    dp_pion =  0.00273402 + (0.00440449)  *1.8    + (-0.00373488)  *pow(1.8,2)    + (0.000996612)  *pow(1.8,3)    + (0) *pow(1.8,4);}
            if(pip_theta > 34 && pip_theta < 35){                     dp_pion =  0.00333542 + (0.00439874)  *pion_p + (-0.00397776)  *pow(pion_p,2) + (0.00105586)   *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 34 && pip_theta < 35 && pion_p >= 1.8){    dp_pion =  0.00333542 + (0.00439874)  *1.8    + (-0.00397776)  *pow(1.8,2)    + (0.00105586)   *pow(1.8,3)    + (0) *pow(1.8,4);}
            if(pip_theta > 35 && pip_theta < 36){                     dp_pion =  0.00354663 + (0.00565397)  *pion_p + (-0.00513503)  *pow(pion_p,2) + (0.00153346)   *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 35 && pip_theta < 36 && pion_p >= 1.8){    dp_pion =  0.00354663 + (0.00565397)  *1.8    + (-0.00513503)  *pow(1.8,2)    + (0.00153346)   *pow(1.8,3)    + (0) *pow(1.8,4);}
            if(pip_theta > 36 && pip_theta < 37){                     dp_pion =  0.00333909 + (0.00842367)  *pion_p + (-0.0077321)   *pow(pion_p,2) + (0.0022489)    *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 36 && pip_theta < 37 && pion_p >= 1.8){    dp_pion =  0.00333909 + (0.00842367)  *1.8    + (-0.0077321)   *pow(1.8,2)    + (0.0022489)    *pow(1.8,3)    + (0) *pow(1.8,4);}
            if(pip_theta > 37 && pip_theta < 38){                     dp_pion =  0.00358828 + (0.0112108)   *pion_p + (-0.0133854)   *pow(pion_p,2) + (0.00486924)   *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 37 && pip_theta < 38 && pion_p >= 1.4){    dp_pion =  0.00358828 + (0.0112108)   *1.4    + (-0.0133854)   *pow(1.4,2)    + (0.00486924)   *pow(1.4,3)    + (0) *pow(1.4,4);}
            if(pip_theta > 38 && pip_theta < 39){                     dp_pion =  0.00354343 + (0.0117121)   *pion_p + (-0.0129649)   *pow(pion_p,2) + (0.00455602)   *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 38 && pip_theta < 39 && pion_p >= 1.3){    dp_pion =  0.00354343 + (0.0117121)   *1.3    + (-0.0129649)   *pow(1.3,2)    + (0.00455602)   *pow(1.3,3)    + (0) *pow(1.3,4);}
            if(pip_theta > 39 && pip_theta < 40){                     dp_pion = -0.00194951 + (0.0409713)   *pion_p + (-0.0595861)   *pow(pion_p,2) + (0.0281588)    *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 39 && pip_theta < 40 && pion_p >= 0.9){    dp_pion = -0.00194951 + (0.0409713)   *0.9    + (-0.0595861)   *pow(0.9,2)    + (0.0281588)    *pow(0.9,3)    + (0) *pow(0.9,4);}
            if(pip_theta > 40 && pip_theta < 41){                     dp_pion = -0.0099217  + (0.0808096)   *pion_p + (-0.119836)    *pow(pion_p,2) + (0.0559553)    *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 40 && pip_theta < 41 && pion_p >= 0.75){   dp_pion = -0.0099217  + (0.0808096)   *0.75   + (-0.119836)    *pow(0.75,2)   + (0.0559553)    *pow(0.75,3)   + (0) *pow(0.75,4);}
            if(pip_theta > 41 && pip_theta < 42){                     dp_pion =  0.00854898 + (0.00025037)  *pion_p + (-0.0113992)   *pow(pion_p,2) + (0.0145178)    *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 41 && pip_theta < 42 && pion_p >= 0.65){   dp_pion =  0.00854898 + (0.00025037)  *0.65   + (-0.0113992)   *pow(0.65,2)   + (0.0145178)    *pow(0.65,3)   + (0) *pow(0.65,4);}
            if(pip_theta > 42){                                       dp_pion =  0.00564818 + (0.00706606)  *pion_p + (0.0042602)    *pow(pion_p,2) + (-0.01141)     *pow(pion_p,3) + (0) *pow(pion_p,4);}
            if(pip_theta > 42 && pion_p >= 0.65){                     dp_pion =  0.00564818 + (0.00706606)  *0.65   + (0.0042602)    *pow(0.65,2)   + (-0.01141)     *pow(0.65,3)   + (0) *pow(0.65,4);}
        }
        if(pion_det == 3){  // Central Detector
            if(pip_theta < 39){                                       dp_pion = -0.045      + (-0.102652) + (0.455589) *pion_p + (-0.671635)   *pow(pion_p,2) + (0.303814)   *pow(pion_p,3);}
            if(pip_theta < 39  && pion_p >= 0.7){                     dp_pion = -0.045      + (-0.102652) + (0.455589) *0.7    + (-0.671635)   *pow(0.7,2)    + (0.303814)   *pow(0.7,3);}
            if(pip_theta > 39  && pip_theta < 40){                    dp_pion =  0.0684552  + (-0.766492)              *pion_p + (1.73092)     *pow(pion_p,2) + (-1.46215)   *pow(pion_p,3) + (0.420127) *pow(pion_p,4);}
            if(pip_theta > 39  && pip_theta < 40 && pion_p >= 1.4){   dp_pion =  0.0684552  + (-0.766492)              *1.4    + (1.73092)     *pow(1.4,2)    + (-1.46215)   *pow(1.4,3)    + (0.420127) *pow(1.4,4);}
            if(pip_theta > 40  && pip_theta < 41){                    dp_pion =  0.751549   + (-7.4593)                *pion_p + (26.8037)     *pow(pion_p,2) + (-47.1576)   *pow(pion_p,3) + (43.8527)  *pow(pion_p,4) + (-20.7039) *pow(pion_p,5) + (3.90931)  *pow(pion_p,6);}
            if(pip_theta > 40  && pip_theta < 41 && pion_p >= 1.45){  dp_pion =  0.751549   + (-7.4593)                *1.45   + (26.8037)     *pow(1.45,2)   + (-47.1576)   *pow(1.45,3)   + (43.8527)  *pow(1.45,4)   + (-20.7039) *pow(1.45,5)   + (3.90931)  *pow(1.45,6);}
            if(pip_theta > 41  && pip_theta < 42){                    dp_pion = -1.35043    + (10.0788)                *pion_p + (-30.4829)    *pow(pion_p,2) + (47.7792)    *pow(pion_p,3) + (-40.996)  *pow(pion_p,4) + (18.2662)  *pow(pion_p,5) + (-3.30449) *pow(pion_p,6);}
            if(pip_theta > 41  && pip_theta < 42 && pion_p >= 1.2){   dp_pion = -1.35043    + (10.0788)                *1.2    + (-30.4829)    *pow(1.2,2)    + (47.7792)    *pow(1.2,3)    + (-40.996)  *pow(1.2,4)    + (18.2662)  *pow(1.2,5)    + (-3.30449) *pow(1.2,6);}
            if(pip_theta > 42  && pip_theta < 43){                    dp_pion = -0.0231195  + (0.0744589)              *pion_p + (-0.0807029)  *pow(pion_p,2) + (0.0264266)  *pow(pion_p,3) + (0)        *pow(pion_p,4);}
            if(pip_theta > 42  && pip_theta < 43 && pion_p >= 1.3){   dp_pion = -0.0231195  + (0.0744589)              *1.3    + (-0.0807029)  *pow(1.3,2)    + (0.0264266)  *pow(1.3,3)    + (0)        *pow(1.3,4);}
            if(pip_theta > 43  && pip_theta < 44){                    dp_pion = -0.00979928 + (0.0351043)              *pion_p + (-0.0365865)  *pow(pion_p,2) + (0.00977218) *pow(pion_p,3) + (0)        *pow(pion_p,4);}
            if(pip_theta > 43  && pip_theta < 44 && pion_p >= 1.1){   dp_pion = -0.00979928 + (0.0351043)              *1.1    + (-0.0365865)  *pow(1.1,2)    + (0.00977218) *pow(1.1,3)    + (0)        *pow(1.1,4);}
            if(pip_theta > 44  && pip_theta < 45){                    dp_pion =  0.00108491 + (-0.00924885)            *pion_p + (0.0216431)   *pow(pion_p,2) + (-0.0137762) *pow(pion_p,3) + (0)        *pow(pion_p,4);}
            if(pip_theta > 44  && pip_theta < 45 && pion_p >= 1.1){   dp_pion =  0.00108491 + (-0.00924885)            *1.1    + (0.0216431)   *pow(1.1,2)    + (-0.0137762) *pow(1.1,3)    + (0)        *pow(1.1,4);}
            if(pip_theta > 45  && pip_theta < 55){                    dp_pion =  0.0092263  + (-0.0676178)             *pion_p + (0.168778)    *pow(pion_p,2) + (-0.167463)  *pow(pion_p,3) + (0.05661)  *pow(pion_p,4);}
            if(pip_theta > 45  && pip_theta < 55 && pion_p >= 1.3){   dp_pion =  0.0092263  + (-0.0676178)             *1.3    + (0.168778)    *pow(1.3,2)    + (-0.167463)  *pow(1.3,3)    + (0.05661)  *pow(1.3,4);}
            if(pip_theta > 55  && pip_theta < 65){                    dp_pion =  0.00805642 + (-0.0670962)             *pion_p + (0.188536)    *pow(pion_p,2) + (-0.20571)   *pow(pion_p,3) + (0.0765)   *pow(pion_p,4);}
            if(pip_theta > 55  && pip_theta < 65 && pion_p >= 1.05){  dp_pion =  0.00805642 + (-0.0670962)             *1.05   + (0.188536)    *pow(1.05,2)   + (-0.20571)   *pow(1.05,3)   + (0.0765)   *pow(1.05,4);}
            if(pip_theta > 65  && pip_theta < 75){                    dp_pion =  0.00312202 + (-0.0269717)             *pion_p + (0.0715236)   *pow(pion_p,2) + (-0.0545622) *pow(pion_p,3) + (0)        *pow(pion_p,4);}
            if(pip_theta > 65  && pip_theta < 75 && pion_p >= 0.75){  dp_pion =  0.00312202 + (-0.0269717)             *0.75   + (0.0715236)   *pow(0.75,2)   + (-0.0545622) *pow(0.75,3)   + (0)        *pow(0.75,4);}
            if(pip_theta > 75  && pip_theta < 85){                    dp_pion =  0.00424971 + (-0.0367683)             *pion_p + (0.10417)     *pow(pion_p,2) + (-0.0899651) *pow(pion_p,3) + (0)        *pow(pion_p,4);}
            if(pip_theta > 75  && pip_theta < 85 && pion_p >= 0.65){  dp_pion =  0.00424971 + (-0.0367683)             *0.65   + (0.10417)     *pow(0.65,2)   + (-0.0899651) *pow(0.65,3)   + (0)        *pow(0.65,4);}
            if(pip_theta > 85  && pip_theta < 95){                    dp_pion =  0.00654123 + (-0.0517915)             *pion_p + (0.147888)    *pow(pion_p,2) + (-0.14253)   *pow(pion_p,3) + (0)        *pow(pion_p,4);}
            if(pip_theta > 85  && pip_theta < 95 && pion_p >= 0.5){   dp_pion =  0.00654123 + (-0.0517915)             *0.5    + (0.147888)    *pow(0.5,2)    + (-0.14253)   *pow(0.5,3)    + (0)        *pow(0.5,4);}
            if(pip_theta > 95  && pip_theta < 105){                   dp_pion = -0.00111721 + (0.00478119)             *pion_p + (0.0158753)   *pow(pion_p,2) + (-0.052902)  *pow(pion_p,3) + (0)        *pow(pion_p,4);}
            if(pip_theta > 95  && pip_theta < 105 && pion_p >= 0.45){ dp_pion = -0.00111721 + (0.00478119)             *0.45   + (0.0158753)   *pow(0.45,2)   + (-0.052902)  *pow(0.45,3)   + (0)        *pow(0.45,4);}
            if(pip_theta > 105 && pip_theta < 115){                   dp_pion = -0.00239839 + (0.00790738)             *pion_p + (0.0311713)   *pow(pion_p,2) + (-0.104157)  *pow(pion_p,3) + (0)        *pow(pion_p,4);}
            if(pip_theta > 105 && pip_theta < 115 && pion_p >= 0.35){ dp_pion = -0.00239839 + (0.00790738)             *0.35   + (0.0311713)   *pow(0.35,2)   + (-0.104157)  *pow(0.35,3)   + (0)        *pow(0.35,4);}
            if(pip_theta > 115 && pip_theta < 125){                   dp_pion = -0.00778793 + (0.0256774)              *pion_p + (0.0932503)   *pow(pion_p,2) + (-0.32771)   *pow(pion_p,3) + (0)        *pow(pion_p,4);}
            if(pip_theta > 115 && pip_theta < 125 && pion_p >= 0.35){ dp_pion = -0.00778793 + (0.0256774)              *0.35   + (0.0932503)   *pow(0.35,2)   + (-0.32771)   *pow(0.35,3)   + (0)        *pow(0.35,4);}
            if(pip_theta > 125 && pip_theta < 135){                   dp_pion = -0.00292778 + (-0.00536697)            *pion_p + (-0.00414351) *pow(pion_p,2) + (0.0196431)  *pow(pion_p,3) + (0)        *pow(pion_p,4);}
            if(pip_theta > 125 && pip_theta < 135 && pion_p >= 0.35){ dp_pion = -0.00292778 + (-0.00536697)            *0.35   + (-0.00414351) *pow(0.35,2)   + (0.0196431)  *pow(0.35,3)   + (0)        *pow(0.35,4);}
        }
    }
    if(outbending == true){ // OUTBENDING
        if(pion_det == 2){  // Forward Detector
            if(pip_theta < 27){                                       dp_pion = 0.00389945  + (-0.004062)    *pion_p + (0.00321842)  *pow(pion_p,2) + (-0.000698299) *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta < 27 && pion_p >= 2.3){                      dp_pion = 0.00389945  + (-0.004062)    *2.3    + (0.00321842)  *pow(2.3,2)    + (-0.000698299) *pow(2.3,3)    + (0)          *pow(2.3,4);}
            if(pip_theta > 27 && pip_theta < 28){                     dp_pion = 0.00727132  + (-0.0117989)   *pion_p + (0.00962999)  *pow(pion_p,2) + (-0.00267005)  *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 27 && pip_theta < 28 && pion_p >= 1.7){    dp_pion = 0.00727132  + (-0.0117989)   *1.7    + (0.00962999)  *pow(1.7,2)    + (-0.00267005)  *pow(1.7,3)    + (0)          *pow(1.7,4);}
            if(pip_theta > 28 && pip_theta < 29){                     dp_pion = 0.00844551  + (-0.0128097)   *pion_p + (0.00945956)  *pow(pion_p,2) + (-0.00237992)  *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 28 && pip_theta < 29 && pion_p >= 2){      dp_pion = 0.00844551  + (-0.0128097)   *2      + (0.00945956)  *pow(2,2)      + (-0.00237992)  *pow(2,3)      + (0)          *pow(2,4);}
            if(pip_theta > 29 && pip_theta < 30){                     dp_pion = 0.00959007  + (-0.0139218)   *pion_p + (0.0122966)   *pow(pion_p,2) + (-0.0034012)   *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 29 && pip_theta < 30 && pion_p >= 1.9){    dp_pion = 0.00959007  + (-0.0139218)   *1.9    + (0.0122966)   *pow(1.9,2)    + (-0.0034012)   *pow(1.9,3)    + (0)          *pow(1.9,4);}
            if(pip_theta > 30 && pip_theta < 31){                     dp_pion = 0.00542816  + (-5.10739e-05) *pion_p + (0.000572038) *pow(pion_p,2) + (-0.000488883) *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 30 && pip_theta < 31 && pion_p >= 1.9){    dp_pion = 0.00542816  + (-5.10739e-05) *1.9    + (0.000572038) *pow(1.9,2)    + (-0.000488883) *pow(1.9,3)    + (0)          *pow(1.9,4);}
            if(pip_theta > 31 && pip_theta < 32){                     dp_pion = 0.0060391   + (-0.000516936) *pion_p + (-0.00286595) *pow(pion_p,2) + (0.00136604)   *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 31 && pip_theta < 32 && pion_p >= 1.8){    dp_pion = 0.0060391   + (-0.000516936) *1.8    + (-0.00286595) *pow(1.8,2)    + (0.00136604)   *pow(1.8,3)    + (0)          *pow(1.8,4);}
            if(pip_theta > 32 && pip_theta < 33){                     dp_pion = 0.0140305   + (-0.0285832)   *pion_p + (0.0248799)   *pow(pion_p,2) + (-0.00701311)  *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 32 && pip_theta < 33 && pion_p >= 1.6){    dp_pion = 0.0140305   + (-0.0285832)   *1.6    + (0.0248799)   *pow(1.6,2)    + (-0.00701311)  *pow(1.6,3)    + (0)          *pow(1.6,4);}
            if(pip_theta > 33 && pip_theta < 34){                     dp_pion = 0.010815    + (-0.0194244)   *pion_p + (0.0174474)   *pow(pion_p,2) + (-0.0049764)   *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 33 && pip_theta < 34 && pion_p >= 1.5){    dp_pion = 0.010815    + (-0.0194244)   *1.5    + (0.0174474)   *pow(1.5,2)    + (-0.0049764)   *pow(1.5,3)    + (0)          *pow(1.5,4);}
            if(pip_theta > 34 && pip_theta < 35){                     dp_pion = 0.0105522   + (-0.0176248)   *pion_p + (0.0161142)   *pow(pion_p,2) + (-0.00472288)  *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 34 && pip_theta < 35 && pion_p >= 1.6){    dp_pion = 0.0105522   + (-0.0176248)   *1.6    + (0.0161142)   *pow(1.6,2)    + (-0.00472288)  *pow(1.6,3)    + (0)          *pow(1.6,4);}
            if(pip_theta > 35 && pip_theta < 36){                     dp_pion = 0.0103938   + (-0.0164003)   *pion_p + (0.0164045)   *pow(pion_p,2) + (-0.00517012)  *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 35 && pip_theta < 36 && pion_p >= 1.5){    dp_pion = 0.0103938   + (-0.0164003)   *1.5    + (0.0164045)   *pow(1.5,2)    + (-0.00517012)  *pow(1.5,3)    + (0)          *pow(1.5,4);}
            if(pip_theta > 36 && pip_theta < 37){                     dp_pion = 0.0441471   + (-0.183937)    *pion_p + (0.338784)    *pow(pion_p,2) + (-0.298985)    *pow(pion_p,3) + (0.126905)   *pow(pion_p,4) + (-0.0208286) *pow(pion_p,5);}
            if(pip_theta > 36 && pip_theta < 37 && pion_p >= 1.8){    dp_pion = 0.0441471   + (-0.183937)    *1.8    + (0.338784)    *pow(1.8,2)    + (-0.298985)    *pow(1.8,3)    + (0.126905)   *pow(1.8,4)    + (-0.0208286) *pow(1.8,5);}
            if(pip_theta > 37 && pip_theta < 38){                     dp_pion = 0.0726119   + (-0.345004)    *pion_p + (0.697789)    *pow(pion_p,2) + (-0.685948)    *pow(pion_p,3) + (0.327195)   *pow(pion_p,4) + (-0.0605621) *pow(pion_p,5);}
            if(pip_theta > 37 && pip_theta < 38 && pion_p >= 1.7){    dp_pion = 0.0726119   + (-0.345004)    *1.7    + (0.697789)    *pow(1.7,2)    + (-0.685948)    *pow(1.7,3)    + (0.327195)   *pow(1.7,4)    + (-0.0605621) *pow(1.7,5);}
            if(pip_theta > 38 && pip_theta < 39){                     dp_pion = 0.0247648   + (-0.0797376)   *pion_p + (0.126535)    *pow(pion_p,2) + (-0.086545)    *pow(pion_p,3) + (0.0219304)  *pow(pion_p,4);}
            if(pip_theta > 38 && pip_theta < 39 && pion_p >= 1.6){    dp_pion = 0.0247648   + (-0.0797376)   *1.6    + (0.126535)    *pow(1.6,2)    + (-0.086545)    *pow(1.6,3)    + (0.0219304)  *pow(1.6,4);}
            if(pip_theta > 39 && pip_theta < 40){                     dp_pion = 0.0208867   + (-0.0492068)   *pion_p + (0.0543187)   *pow(pion_p,2) + (-0.0183393)   *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 39 && pip_theta < 40 && pion_p >= 1.2){    dp_pion = 0.0208867   + (-0.0492068)   *1.2    + (0.0543187)   *pow(1.2,2)    + (-0.0183393)   *pow(1.2,3)    + (0)          *pow(1.2,4);}
            if(pip_theta > 40 && pip_theta < 41){                     dp_pion = 0.0148655   + (-0.0203483)   *pion_p + (0.00835867)  *pow(pion_p,2) + (0.00697134)   *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 40 && pip_theta < 41 && pion_p >= 1.0){    dp_pion = 0.0148655   + (-0.0203483)   *1.0    + (0.00835867)  *pow(1.0,2)    + (0.00697134)   *pow(1.0,3)    + (0)          *pow(1.0,4);}
            if(pip_theta > 41 && pip_theta < 42){                     dp_pion = 0.0223585   + (-0.0365262)   *pion_p + (-0.0150027)  *pow(pion_p,2) + (0.0854164)    *pow(pion_p,3) + (-0.0462718) *pow(pion_p,4);}
            if(pip_theta > 41 && pip_theta < 42 && pion_p >= 0.7){    dp_pion = 0.007617;}
            if(pip_theta > 42){                                       dp_pion = 0.0152373   + (-0.0106377)   *pion_p + (-0.0257573)  *pow(pion_p,2) + (0.0344851)    *pow(pion_p,3) + (0)          *pow(pion_p,4);}
            if(pip_theta > 42 && pion_p >= 0.75){                     dp_pion = 0.0152373   + (-0.0106377)   *0.75   + (-0.0257573)  *pow(0.75,2)   + (0.0344851)    *pow(0.75,3)   + (0)          *pow(0.75,4);}
        }
        if(pion_det == 3){ // Central Detector
            if(pip_theta < 39){                                       dp_pion = -0.05        + (-0.0758897) + (0.362231) *pion_p + (-0.542404)   *pow(pion_p,2) + (0.241344)   *pow(pion_p,3);}
            if(pip_theta < 39  && pion_p >= 0.8){                     dp_pion = -0.05        + (-0.0758897) + (0.362231) *0.8    + (-0.542404)   *pow(0.8,2)    + (0.241344)   *pow(0.8,3);}
            if(pip_theta > 39  && pip_theta < 40){                    dp_pion =  0.0355259   + (-0.589712)               *pion_p + (1.4206)      *pow(pion_p,2) + (-1.24179)   *pow(pion_p,3) + (0.365524)  *pow(pion_p,4);}
            if(pip_theta > 39  && pip_theta < 40  && pion_p >= 1.35){ dp_pion =  0.0355259   + (-0.589712)               *1.35   + (1.4206)      *pow(1.35,2)   + (-1.24179)   *pow(1.35,3)   + (0.365524)  *pow(1.35,4);}
            if(pip_theta > 40  && pip_theta < 41){                    dp_pion = -0.252336    + (1.02032)                 *pion_p + (-1.51461)    *pow(pion_p,2) + (0.967772)   *pow(pion_p,3) + (-0.226028) *pow(pion_p,4);}
            if(pip_theta > 40  && pip_theta < 41  && pion_p >= 1.4){  dp_pion = -0.252336    + (1.02032)                 *1.4    + (-1.51461)    *pow(1.4,2)    + (0.967772)   *pow(1.4,3)    + (-0.226028) *pow(1.4,4);}
            if(pip_theta > 41  && pip_theta < 42){                    dp_pion = -0.710129    + (4.49613)                 *pion_p + (-11.01)      *pow(pion_p,2) + (12.9945)    *pow(pion_p,3) + (-7.41641)  *pow(pion_p,4) + (1.63923)   *pow(pion_p,5);}
            if(pip_theta > 41  && pip_theta < 42  && pion_p >= 1.2){  dp_pion = -0.710129    + (4.49613)                 *1.2    + (-11.01)      *pow(1.2,2)    + (12.9945)    *pow(1.2,3)    + (-7.41641)  *pow(1.2,4)    + (1.63923)   *pow(1.2,5);}
            if(pip_theta > 42  && pip_theta < 43){                    dp_pion = -0.0254912   + (0.0851432)               *pion_p + (-0.0968583)  *pow(pion_p,2) + (0.0350334)  *pow(pion_p,3) + (0)         *pow(pion_p,4);}
            if(pip_theta > 42  && pip_theta < 43  && pion_p >= 1.2){  dp_pion = -0.0254912   + (0.0851432)               *1.2    + (-0.0968583)  *pow(1.2,2)    + (0.0350334)  *pow(1.2,3)    + (0)         *pow(1.2,4);}
            if(pip_theta > 43  && pip_theta < 44){                    dp_pion = -0.0115965   + (0.0438726)               *pion_p + (-0.0500474)  *pow(pion_p,2) + (0.0163627)  *pow(pion_p,3) + (0)         *pow(pion_p,4);}
            if(pip_theta > 43  && pip_theta < 44  && pion_p >= 1.4){  dp_pion = -0.0115965   + (0.0438726)               *1.4    + (-0.0500474)  *pow(1.4,2)    + (0.0163627)  *pow(1.4,3)    + (0)         *pow(1.4,4);}
            if(pip_theta > 44  && pip_theta < 45){                    dp_pion =  0.00273414  + (-0.01851)                *pion_p + (0.0377032)   *pow(pion_p,2) + (-0.0226696) *pow(pion_p,3) + (0)         *pow(pion_p,4);}
            if(pip_theta > 44  && pip_theta < 45  && pion_p >= 1){    dp_pion =  0.00273414  + (-0.01851)                *1      + (0.0377032)   *pow(1,2)      + (-0.0226696) *pow(1,3)      + (0)         *pow(1,4);}
            if(pip_theta > 45  && pip_theta < 55){                    dp_pion =  0.0271952   + (-0.25981)                *pion_p + (0.960051)    *pow(pion_p,2) + (-1.76651)   *pow(pion_p,3) + (1.72872)   *pow(pion_p,4) + (-0.856946) *pow(pion_p,5) + (0.167564) *pow(pion_p,6);}
            if(pip_theta > 45  && pip_theta < 55  && pion_p >= 1.4){  dp_pion =  0.0271952   + (-0.25981)                *1.4    + (0.960051)    *pow(1.4,2)    + (-1.76651)   *pow(1.4,3)    + (1.72872)   *pow(1.4,4)    + (-0.856946) *pow(1.4,5)    + (0.167564) *pow(1.4,6);}
            if(pip_theta > 55  && pip_theta < 65){                    dp_pion =  0.00734975  + (-0.0598841)              *pion_p + (0.161495)    *pow(pion_p,2) + (-0.1629)    *pow(pion_p,3) + (0.0530098) *pow(pion_p,4);}
            if(pip_theta > 55  && pip_theta < 65  && pion_p >= 1.2){  dp_pion =  0.00734975  + (-0.0598841)              *1.2    + (0.161495)    *pow(1.2,2)    + (-0.1629)    *pow(1.2,3)    + (0.0530098) *pow(1.2,4);}
            if(pip_theta > 65  && pip_theta < 75){                    dp_pion =  0.00321351  + (-0.0289322)              *pion_p + (0.0786484)   *pow(pion_p,2) + (-0.0607041) *pow(pion_p,3) + (0)         *pow(pion_p,4);}
            if(pip_theta > 65  && pip_theta < 75  && pion_p >= 0.95){ dp_pion =  0.00321351  + (-0.0289322)              *0.95   + (0.0786484)   *pow(0.95,2)   + (-0.0607041) *pow(0.95,3)   + (0)         *pow(0.95,4);}
            if(pip_theta > 75  && pip_theta < 85){                    dp_pion =  0.00644253  + (-0.0543896)              *pion_p + (0.148933)    *pow(pion_p,2) + (-0.1256)    *pow(pion_p,3) + (0)         *pow(pion_p,4);}
            if(pip_theta > 75  && pip_theta < 85  && pion_p >= 0.7){  dp_pion =  0.00644253  + (-0.0543896)              *0.7    + (0.148933)    *pow(0.7,2)    + (-0.1256)    *pow(0.7,3)    + (0)         *pow(0.7,4);}
            if(pip_theta > 85  && pip_theta < 95){                    dp_pion =  0.00671152  + (-0.0537269)              *pion_p + (0.154509)    *pow(pion_p,2) + (-0.147667)  *pow(pion_p,3) + (0)         *pow(pion_p,4);}
            if(pip_theta > 85  && pip_theta < 95  && pion_p >= 0.65){ dp_pion =  0.00671152  + (-0.0537269)              *0.65   + (0.154509)    *pow(0.65,2)   + (-0.147667)  *pow(0.65,3)   + (0)         *pow(0.65,4);}
            if(pip_theta > 95  && pip_theta < 105){                   dp_pion = -0.000709077 + (0.00331818)              *pion_p + (0.0109241)   *pow(pion_p,2) + (-0.0351682) *pow(pion_p,3) + (0)         *pow(pion_p,4);}
            if(pip_theta > 95  && pip_theta < 105 && pion_p >= 0.45){ dp_pion = -0.000709077 + (0.00331818)              *0.45   + (0.0109241)   *pow(0.45,2)   + (-0.0351682) *pow(0.45,3)   + (0)         *pow(0.45,4);}
            if(pip_theta > 105 && pip_theta < 115){                   dp_pion = -0.00260164  + (0.00846919)              *pion_p + (0.0315497)   *pow(pion_p,2) + (-0.105756)  *pow(pion_p,3) + (0)         *pow(pion_p,4);}
            if(pip_theta > 105 && pip_theta < 115 && pion_p >= 0.45){ dp_pion = -0.00260164  + (0.00846919)              *0.45   + (0.0315497)   *pow(0.45,2)   + (-0.105756)  *pow(0.45,3)   + (0)         *pow(0.45,4);}
            if(pip_theta > 115 && pip_theta < 125){                   dp_pion = -0.00544336  + (0.018256)                *pion_p + (0.0664618)   *pow(pion_p,2) + (-0.240312)  *pow(pion_p,3) + (0)         *pow(pion_p,4);}
            if(pip_theta > 115 && pip_theta < 125 && pion_p >= 0.45){ dp_pion = -0.00544336  + (0.018256)                *0.45   + (0.0664618)   *pow(0.45,2)   + (-0.240312)  *pow(0.45,3)   + (0)         *pow(0.45,4);}
            if(pip_theta > 125 && pip_theta < 135){                   dp_pion = -0.00281073  + (-0.00495863)             *pion_p + (-0.00362356) *pow(pion_p,2) + (0.0178764)  *pow(pion_p,3) + (0)         *pow(pion_p,4);}
            if(pip_theta > 125 && pip_theta < 135 && pion_p >= 0.35){ dp_pion = -0.00281073  + (-0.00495863)             *0.35   + (-0.00362356) *pow(0.35,2)   + (0.0178764)  *pow(0.35,3)   + (0)         *pow(0.35,4);}
        }
    }

    return dp_pion;
};

