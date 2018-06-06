// Copyright (C) 2011-2017 Vincent Heuveline
//
// HiFlow3 is free software: you can redistribute it and/or modify it under the
// terms of the European Union Public Licence (EUPL) v1.2 as published by the
// European Union or (at your option) any later version.
//
// HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for more
// details.
//
// You should have received a copy of the European Union Public Licence (EUPL) v1.2
// along with HiFlow3.  If not, see <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

// $Id: lagrangecoeff11.cc,v 1.1 2000/03/11 14:14:07 heuvelin Exp $

/// \author Martin Baumann, Michael Schick

#include "lagrangecoeff.h"
#include "lagrangecoeff11.h"

namespace hiflow
{

    static double _poly_lagrange_coeff11[] = {
                                              1.000000000000000000000000000000e+00,
                                              -3.321865079365079365079365079370e+01,
                                              4.574784325396825396825396825400e+02,
                                              -3.509762785493827160493827160490e+03,
                                              1.687067808366402116402116402120e+04,
                                              -5.382058309496252204585537918870e+04,
                                              1.170583535763888888888888888890e+05,
                                              -1.744920214128637566137566137570e+05,
                                              1.754425365327380952380952380950e+05,
                                              -1.137127551601080246913580246910e+05,
                                              4.288595337466931216931216931220e+04,
                                              -7.147658895778218694885361552030e+03,

                                              0.000000000000000000000000000000e-01,
                                              1.210000000000000000000000000000e+02,
                                              -2.688456746031746031746031746030e+03,
                                              2.578186613095238095238095238100e+04,
                                              -1.410807696042768959435626102290e+05,
                                              4.894635824763007054673721340390e+05,
                                              -1.128191147251157407407407407410e+06,
                                              1.753958162980324074074074074070e+06,
                                              -1.819994798172949735449735449740e+06,
                                              1.208604140558862433862433862430e+06,
                                              -4.645978282255842151675485008820e+05,
                                              7.862424785356040564373897707230e+04,

                                              0.000000000000000000000000000000e-01,
                                              -3.025000000000000000000000000000e+02,
                                              8.384891865079365079365079365080e+03,
                                              -9.227032058531746031746031746030e+04,
                                              5.542164793926366843033509700180e+05,
                                              -2.055189483648864638447971781310e+06,
                                              4.977184226157407407407407407410e+06,
                                              -8.035638712991898148148148148150e+06,
                                              8.587823555935846560846560846560e+06,
                                              -5.838337743506117724867724867720e+06,
                                              2.287250846649029982363315696650e+06,
                                              -3.931212392678020282186948853620e+05,

                                              0.000000000000000000000000000000e-01,
                                              6.050000000000000000000000000000e+02,
                                              -1.787895039682539682539682539680e+04,
                                              2.112183002314814814814814814810e+05,
                                              -1.348939384375000000000000000000e+06,
                                              5.260649164575066137566137566140e+06,
                                              -1.327240583567708333333333333330e+07,
                                              2.215481584956597222222222222220e+07,
                                              -2.433334817304067460317460317460e+07,
                                              1.692045796782407407407407407410e+07,
                                              -6.754537656510416666666666666670e+06,
                                              1.179363717803406084656084656080e+06,

                                              0.000000000000000000000000000000e-01,
                                              -9.075000000000000000000000000000e+02,
                                              2.765030059523809523809523809520e+04,
                                              -3.391233508928571428571428571430e+05,
                                              2.252520512880291005291005291010e+06,
                                              -9.115708950504298941798941798940e+06,
                                              2.377397954479166666666666666670e+07,
                                              -4.085201212239583333333333333330e+07,
                                              4.600847609558531746031746031750e+07,
                                              -3.269079264060019841269841269840e+07,
                                              1.329464554614748677248677248680e+07,
                                              -2.358727435606812169312169312170e+06,

                                              0.000000000000000000000000000000e-01,
                                              1.016400000000000000000000000000e+03,
                                              -3.152735666666666666666666666670e+04,
                                              3.956208941666666666666666666670e+05,
                                              -2.696956928009259259259259259260e+06,
                                              1.121405196261574074074074074070e+07,
                                              -3.003232633996527777777777777780e+07,
                                              5.290699262711805555555555555560e+07,
                                              -6.095830678437500000000000000000e+07,
                                              4.421151920625000000000000000000e+07,
                                              -1.831230209098379629629629629630e+07,
                                              3.302218409849537037037037037040e+06,

                                              0.000000000000000000000000000000e-01,
                                              -8.470000000000000000000000000000e+02,
                                              2.658336388888888888888888888890e+04,
                                              -3.387480652314814814814814814810e+05,
                                              2.351730959722222222222222222220e+06,
                                              -9.977957577372685185185185185190e+06,
                                              2.729311165625000000000000000000e+07,
                                              -4.911105410940972222222222222220e+07,
                                              5.775780960277777777777777777780e+07,
                                              -4.271051083813657407407407407410e+07,
                                              1.801210041736111111111111111110e+07,
                                              -3.302218409849537037037037037040e+06,

                                              0.000000000000000000000000000000e-01,
                                              5.185714285714285714285714285710e+02,
                                              -1.641134523809523809523809523810e+04,
                                              2.114459875000000000000000000000e+05,
                                              -1.487790435548941798941798941800e+06,
                                              6.410695236094576719576719576720e+06,
                                              -1.783586700538194444444444444440e+07,
                                              3.267532663187003968253968253970e+07,
                                              -3.913963496830357142857142857140e+07,
                                              2.947434613750000000000000000000e+07,
                                              -1.265135624552744708994708994710e+07,
                                              2.358727435606812169312169312170e+06,

                                              0.000000000000000000000000000000e-01,
                                              -2.268750000000000000000000000000e+02,
                                              7.224528273809523809523809523810e+03,
                                              -9.385669300595238095238095238100e+04,
                                              6.672244790757275132275132275130e+05,
                                              -2.910101431502149470899470899470e+06,
                                              8.209155321354166666666666666670e+06,
                                              -1.527002540078125000000000000000e+07,
                                              1.859159243196924603174603174600e+07,
                                              -1.424008588190724206349206349210e+07,
                                              6.218463239327050264550264550260e+06,
                                              -1.179363717803406084656084656080e+06,

                                              0.000000000000000000000000000000e-01,
                                              6.722222222222222222222222222220e+01,
                                              -2.150871031746031746031746031750e+03,
                                              2.812387448192239858906525573190e+04,
                                              -2.015604295469576719576719576720e+05,
                                              8.877328350666887125220458553790e+05,
                                              -2.532932398524305555555555555560e+06,
                                              4.773116392216435185185185185190e+06,
                                              -5.895932515600198412698412698410e+06,
                                              4.587497436744929453262786596120e+06,
                                              -2.037082785296792328042328042330e+06,
                                              3.931212392678020282186948853620e+05,

                                              0.000000000000000000000000000000e-01,
                                              -1.210000000000000000000000000000e+01,
                                              3.886356746031746031746031746030e+02,
                                              -5.107989791666666666666666666670e+03,
                                              3.684934093364197530864197530860e+04,
                                              -1.636009297853284832451499118170e+05,
                                              4.712680326851851851851851851850e+05,
                                              -8.980112423206018518518518518520e+05,
                                              1.123541092542989417989417989420e+06,
                                              -8.869594902488425925925925925930e+05,
                                              4.002688981635802469135802469140e+05,
                                              -7.862424785356040564373897707230e+04,

                                              0.000000000000000000000000000000e-01,
                                              1.000000000000000000000000000000e+00,
                                              -3.221865079365079365079365079370e+01,
                                              4.252597817460317460317460317460e+02,
                                              -3.084503003747795414462081128750e+03,
                                              1.378617507991622574955908289240e+04,
                                              -4.003440801504629629629629629630e+04,
                                              7.702394556134259259259259259260e+04,
                                              -9.746807585152116402116402116400e+04,
                                              7.797446068121693121693121693120e+04,
                                              -3.573829447889109347442680776010e+04,
                                              7.147658895778218694885361552030e+03

    };

    static double _poly_x_lagrange_coeff11[] = {
                                                -3.321865079365079365079365079370e+01,
                                                9.149568650793650793650793650790e+02,
                                                -1.052928835648148148148148148150e+04,
                                                6.748271233465608465608465608470e+04,
                                                -2.691029154748126102292768959440e+05,
                                                7.023501214583333333333333333330e+05,
                                                -1.221444149890046296296296296300e+06,
                                                1.403540292261904761904761904760e+06,
                                                -1.023414796440972222222222222220e+06,
                                                4.288595337466931216931216931220e+05,
                                                -7.862424785356040564373897707230e+04,
                                                0.000000000000000000000000000000e-01,

                                                1.210000000000000000000000000000e+02,
                                                -5.376913492063492063492063492060e+03,
                                                7.734559839285714285714285714290e+04,
                                                -5.643230784171075837742504409170e+05,
                                                2.447317912381503527336860670190e+06,
                                                -6.769146883506944444444444444440e+06,
                                                1.227770714086226851851851851850e+07,
                                                -1.455995838538359788359788359790e+07,
                                                1.087743726502976190476190476190e+07,
                                                -4.645978282255842151675485008820e+06,
                                                8.648667263891644620811287477950e+05,
                                                0.000000000000000000000000000000e-01,

                                                -3.025000000000000000000000000000e+02,
                                                1.676978373015873015873015873020e+04,
                                                -2.768109617559523809523809523810e+05,
                                                2.216865917570546737213403880070e+06,
                                                -1.027594741824432319223985890650e+07,
                                                2.986310535694444444444444444440e+07,
                                                -5.624947099094328703703703703700e+07,
                                                6.870258844748677248677248677250e+07,
                                                -5.254503969155505952380952380950e+07,
                                                2.287250846649029982363315696650e+07,
                                                -4.324333631945822310405643738980e+06,
                                                0.000000000000000000000000000000e-01,

                                                6.050000000000000000000000000000e+02,
                                                -3.575790079365079365079365079370e+04,
                                                6.336549006944444444444444444440e+05,
                                                -5.395757537500000000000000000000e+06,
                                                2.630324582287533068783068783070e+07,
                                                -7.963443501406250000000000000000e+07,
                                                1.550837109469618055555555555560e+08,
                                                -1.946667853843253968253968253970e+08,
                                                1.522841217104166666666666666670e+08,
                                                -6.754537656510416666666666666670e+07,
                                                1.297300089583746693121693121690e+07,
                                                0.000000000000000000000000000000e-01,

                                                -9.075000000000000000000000000000e+02,
                                                5.530060119047619047619047619050e+04,
                                                -1.017370052678571428571428571430e+06,
                                                9.010082051521164021164021164020e+06,
                                                -4.557854475252149470899470899470e+07,
                                                1.426438772687500000000000000000e+08,
                                                -2.859640848567708333333333333330e+08,
                                                3.680678087646825396825396825400e+08,
                                                -2.942171337654017857142857142860e+08,
                                                1.329464554614748677248677248680e+08,
                                                -2.594600179167493386243386243390e+07,
                                                0.000000000000000000000000000000e-01,

                                                1.016400000000000000000000000000e+03,
                                                -6.305471333333333333333333333330e+04,
                                                1.186862682500000000000000000000e+06,
                                                -1.078782771203703703703703703700e+07,
                                                5.607025981307870370370370370370e+07,
                                                -1.801939580397916666666666666670e+08,
                                                3.703489483898263888888888888890e+08,
                                                -4.876664542750000000000000000000e+08,
                                                3.979036728562500000000000000000e+08,
                                                -1.831230209098379629629629629630e+08,
                                                3.632440250834490740740740740740e+07,
                                                0.000000000000000000000000000000e-01,

                                                -8.470000000000000000000000000000e+02,
                                                5.316672777777777777777777777780e+04,
                                                -1.016244195694444444444444444440e+06,
                                                9.406923838888888888888888888890e+06,
                                                -4.988978788686342592592592592590e+07,
                                                1.637586699375000000000000000000e+08,
                                                -3.437773787658680555555555555560e+08,
                                                4.620624768222222222222222222220e+08,
                                                -3.843945975432291666666666666670e+08,
                                                1.801210041736111111111111111110e+08,
                                                -3.632440250834490740740740740740e+07,
                                                0.000000000000000000000000000000e-01,

                                                5.185714285714285714285714285710e+02,
                                                -3.282269047619047619047619047620e+04,
                                                6.343379625000000000000000000000e+05,
                                                -5.951161742195767195767195767200e+06,
                                                3.205347618047288359788359788360e+07,
                                                -1.070152020322916666666666666670e+08,
                                                2.287272864230902777777777777780e+08,
                                                -3.131170797464285714285714285710e+08,
                                                2.652691152375000000000000000000e+08,
                                                -1.265135624552744708994708994710e+08,
                                                2.594600179167493386243386243390e+07,
                                                0.000000000000000000000000000000e-01,

                                                -2.268750000000000000000000000000e+02,
                                                1.444905654761904761904761904760e+04,
                                                -2.815700790178571428571428571430e+05,
                                                2.668897916302910052910052910050e+06,
                                                -1.455050715751074735449735449740e+07,
                                                4.925493192812500000000000000000e+07,
                                                -1.068901778054687500000000000000e+08,
                                                1.487327394557539682539682539680e+08,
                                                -1.281607729371651785714285714290e+08,
                                                6.218463239327050264550264550260e+07,
                                                -1.297300089583746693121693121690e+07,
                                                0.000000000000000000000000000000e-01,

                                                6.722222222222222222222222222220e+01,
                                                -4.301742063492063492063492063490e+03,
                                                8.437162344576719576719576719580e+04,
                                                -8.062417181878306878306878306880e+05,
                                                4.438664175333443562610229276900e+06,
                                                -1.519759439114583333333333333330e+07,
                                                3.341181474551504629629629629630e+07,
                                                -4.716746012480158730158730158730e+07,
                                                4.128747693070436507936507936510e+07,
                                                -2.037082785296792328042328042330e+07,
                                                4.324333631945822310405643738980e+06,
                                                0.000000000000000000000000000000e-01,

                                                -1.210000000000000000000000000000e+01,
                                                7.772713492063492063492063492060e+02,
                                                -1.532396937500000000000000000000e+04,
                                                1.473973637345679012345679012350e+05,
                                                -8.180046489266424162257495590830e+05,
                                                2.827608196111111111111111111110e+06,
                                                -6.286078696244212962962962962960e+06,
                                                8.988328740343915343915343915340e+06,
                                                -7.982635412239583333333333333330e+06,
                                                4.002688981635802469135802469140e+06,
                                                -8.648667263891644620811287477950e+05,
                                                0.000000000000000000000000000000e-01,

                                                1.000000000000000000000000000000e+00,
                                                -6.443730158730158730158730158730e+01,
                                                1.275779345238095238095238095240e+03,
                                                -1.233801201499118165784832451500e+04,
                                                6.893087539958112874779541446210e+04,
                                                -2.402064480902777777777777777780e+05,
                                                5.391676189293981481481481481480e+05,
                                                -7.797446068121693121693121693120e+05,
                                                7.017701461309523809523809523810e+05,
                                                -3.573829447889109347442680776010e+05,
                                                7.862424785356040564373897707230e+04,
                                                0.000000000000000000000000000000e-01

    };

    static double _poly_xx_lagrange_coeff11[] = {
                                                 9.149568650793650793650793650790e+02,
                                                 -2.105857671296296296296296296300e+04,
                                                 2.024481370039682539682539682540e+05,
                                                 -1.076411661899250440917107583770e+06,
                                                 3.511750607291666666666666666670e+06,
                                                 -7.328664899340277777777777777780e+06,
                                                 9.824782045833333333333333333330e+06,
                                                 -8.187318371527777777777777777780e+06,
                                                 3.859735803720238095238095238100e+06,
                                                 -7.862424785356040564373897707230e+05,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01,

                                                 -5.376913492063492063492063492060e+03,
                                                 1.546911967857142857142857142860e+05,
                                                 -1.692969235251322751322751322750e+06,
                                                 9.789271649526014109347442680780e+06,
                                                 -3.384573441753472222222222222220e+07,
                                                 7.366624284517361111111111111110e+07,
                                                 -1.019197086976851851851851851850e+08,
                                                 8.701949812023809523809523809520e+07,
                                                 -4.181380454030257936507936507940e+07,
                                                 8.648667263891644620811287477950e+06,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01,

                                                 1.676978373015873015873015873020e+04,
                                                 -5.536219235119047619047619047620e+05,
                                                 6.650597752711640211640211640210e+06,
                                                 -4.110378967297729276895943562610e+07,
                                                 1.493155267847222222222222222220e+08,
                                                 -3.374968259456597222222222222220e+08,
                                                 4.809181191324074074074074074070e+08,
                                                 -4.203603175324404761904761904760e+08,
                                                 2.058525761984126984126984126980e+08,
                                                 -4.324333631945822310405643738980e+07,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01,

                                                 -3.575790079365079365079365079370e+04,
                                                 1.267309801388888888888888888890e+06,
                                                 -1.618727261250000000000000000000e+07,
                                                 1.052129832915013227513227513230e+08,
                                                 -3.981721750703125000000000000000e+08,
                                                 9.305022656817708333333333333330e+08,
                                                 -1.362667497690277777777777777780e+09,
                                                 1.218272973683333333333333333330e+09,
                                                 -6.079083890859375000000000000000e+08,
                                                 1.297300089583746693121693121690e+08,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01,

                                                 5.530060119047619047619047619050e+04,
                                                 -2.034740105357142857142857142860e+06,
                                                 2.703024615456349206349206349210e+07,
                                                 -1.823141790100859788359788359790e+08,
                                                 7.132193863437500000000000000000e+08,
                                                 -1.715784509140625000000000000000e+09,
                                                 2.576474661352777777777777777780e+09,
                                                 -2.353737070123214285714285714290e+09,
                                                 1.196518099153273809523809523810e+09,
                                                 -2.594600179167493386243386243390e+08,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01,

                                                 -6.305471333333333333333333333330e+04,
                                                 2.373725365000000000000000000000e+06,
                                                 -3.236348313611111111111111111110e+07,
                                                 2.242810392523148148148148148150e+08,
                                                 -9.009697901989583333333333333330e+08,
                                                 2.222093690338958333333333333330e+09,
                                                 -3.413665179925000000000000000000e+09,
                                                 3.183229382850000000000000000000e+09,
                                                 -1.648107188188541666666666666670e+09,
                                                 3.632440250834490740740740740740e+08,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01,

                                                 5.316672777777777777777777777780e+04,
                                                 -2.032488391388888888888888888890e+06,
                                                 2.822077151666666666666666666670e+07,
                                                 -1.995591515474537037037037037040e+08,
                                                 8.187933496875000000000000000000e+08,
                                                 -2.062664272595208333333333333330e+09,
                                                 3.234437337755555555555555555560e+09,
                                                 -3.075156780345833333333333333330e+09,
                                                 1.621089037562500000000000000000e+09,
                                                 -3.632440250834490740740740740740e+08,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01,

                                                 -3.282269047619047619047619047620e+04,
                                                 1.268675925000000000000000000000e+06,
                                                 -1.785348522658730158730158730160e+07,
                                                 1.282139047218915343915343915340e+08,
                                                 -5.350760101614583333333333333330e+08,
                                                 1.372363718538541666666666666670e+09,
                                                 -2.191819558225000000000000000000e+09,
                                                 2.122152921900000000000000000000e+09,
                                                 -1.138622062097470238095238095240e+09,
                                                 2.594600179167493386243386243390e+08,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01,

                                                 1.444905654761904761904761904760e+04,
                                                 -5.631401580357142857142857142860e+05,
                                                 8.006693748908730158730158730160e+06,
                                                 -5.820202863004298941798941798940e+07,
                                                 2.462746596406250000000000000000e+08,
                                                 -6.413410668328125000000000000000e+08,
                                                 1.041129176190277777777777777780e+09,
                                                 -1.025286183497321428571428571430e+09,
                                                 5.596616915394345238095238095240e+08,
                                                 -1.297300089583746693121693121690e+08,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01,

                                                 -4.301742063492063492063492063490e+03,
                                                 1.687432468915343915343915343920e+05,
                                                 -2.418725154563492063492063492060e+06,
                                                 1.775465670133377425044091710760e+07,
                                                 -7.598797195572916666666666666670e+07,
                                                 2.004708884730902777777777777780e+08,
                                                 -3.301722208736111111111111111110e+08,
                                                 3.302998154456349206349206349210e+08,
                                                 -1.833374506767113095238095238100e+08,
                                                 4.324333631945822310405643738980e+07,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01,

                                                 7.772713492063492063492063492060e+02,
                                                 -3.064793875000000000000000000000e+04,
                                                 4.421920912037037037037037037040e+05,
                                                 -3.272018595706569664902998236330e+06,
                                                 1.413804098055555555555555555560e+07,
                                                 -3.771647217746527777777777777780e+07,
                                                 6.291830118240740740740740740740e+07,
                                                 -6.386108329791666666666666666670e+07,
                                                 3.602420083472222222222222222220e+07,
                                                 -8.648667263891644620811287477950e+06,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01,

                                                 -6.443730158730158730158730158730e+01,
                                                 2.551558690476190476190476190480e+03,
                                                 -3.701403604497354497354497354500e+04,
                                                 2.757235015983245149911816578480e+05,
                                                 -1.201032240451388888888888888890e+06,
                                                 3.235005713576388888888888888890e+06,
                                                 -5.458212247685185185185185185190e+06,
                                                 5.614161169047619047619047619050e+06,
                                                 -3.216446503100198412698412698410e+06,
                                                 7.862424785356040564373897707230e+05,
                                                 0.000000000000000000000000000000e-01,
                                                 0.000000000000000000000000000000e-01
    };

    LagrangeCoeff _lagrange_coeff11 ( 11,
                                      _poly_lagrange_coeff11,
                                      _poly_x_lagrange_coeff11,
                                      _poly_xx_lagrange_coeff11 );
} // namespace hiflow
