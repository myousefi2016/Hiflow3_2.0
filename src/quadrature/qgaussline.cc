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

/// \author Staffan Ronnas

#include "qgaussline.h"

namespace hiflow
{

    template<class DataType>
    QuadratureGaussLine<DataType>::QuadratureGaussLine ( ) : QuadratureType<DataType>::QuadratureType ( )
    {
        //////////////////////////////////////////////////////////////
        // Implemented sizes: (please fill up if a new size was added)
        // 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
        //////////////////////////////////////////////////////////////
        int total_number_of_quadratures = 15;

        int cntr = 0;
        int size = 0;

        // Create order->size map.
        order_size_map_.clear ( );
        order_size_map_.resize ( 30 );
        order_size_map_[0] = 1;
        order_size_map_[1] = 1;
        order_size_map_[2] = 2;
        order_size_map_[3] = 2;
        order_size_map_[4] = 3;
        order_size_map_[5] = 3;
        order_size_map_[6] = 4;
        order_size_map_[7] = 4;
        order_size_map_[8] = 5;
        order_size_map_[9] = 5;

        order_size_map_[10] = 6;
        order_size_map_[11] = 6;
        order_size_map_[12] = 7;
        order_size_map_[13] = 7;
        order_size_map_[14] = 8;
        order_size_map_[15] = 8;
        order_size_map_[16] = 9;
        order_size_map_[17] = 9;
        order_size_map_[18] = 10;
        order_size_map_[19] = 10;

        order_size_map_[20] = 11;
        order_size_map_[21] = 11;
        order_size_map_[22] = 12;
        order_size_map_[23] = 12;
        order_size_map_[24] = 13;
        order_size_map_[25] = 13;
        order_size_map_[26] = 14;
        order_size_map_[27] = 14;
        order_size_map_[28] = 15;
        order_size_map_[29] = 15;

        // First resize vector fields to the number of implemented quadratures
        x_.resize ( total_number_of_quadratures );
        y_.resize ( total_number_of_quadratures );
        z_.resize ( total_number_of_quadratures );
        w_.resize ( total_number_of_quadratures );

        // Next fill vector fields with quadrature data

        // ---------------------------------------------
        size = 1;

        DataType x1[] = {
                         .5
        };
        DataType w1[] = {
                         1.0
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x1, x1 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w1, w1 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 2;

        DataType x2[] = {
                         .211324865405187117745425609748,
                         .788675134594812882254574390252
        };
        DataType w2[] = {
                         .5,
                         .5
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x2, x2 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w2, w2 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 3;

        DataType x3[] = {
                         .112701665379258311482073460022,
                         .5,
                         .887298334620741688517926539978
        };
        DataType w3[] = {
                         .277777777777777777777777777779,
                         .444444444444444444444444444445,
                         .277777777777777777777777777779
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x3, x3 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w3, w3 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 4;

        DataType x4[] = {
                         .069431844202973712388026755553,
                         .330009478207571867598667120448,
                         .669990521792428132401332879552,
                         .930568155797026287611973244447
        };
        DataType w4[] = {
                         .173927422568726928686531974611,
                         .326072577431273071313468025388,
                         .326072577431273071313468025388,
                         .173927422568726928686531974611
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x4, x4 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w4, w4 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 5;

        DataType x5[] = {
                         .046910077030668003601186560850,
                         .23076534494715845448184278965,
                         .5,
                         .76923465505284154551815721035,
                         .95308992296933199639881343915
        };
        DataType w5[] = {
                         .118463442528094543757132020362,
                         .239314335249683234020645757416,
                         .284444444444444444444444444445,
                         .239314335249683234020645757416,
                         .118463442528094543757132020362
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x5, x5 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w5, w5 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 6;

        DataType x6[] = {
                         .033765242898423986093849222756,
                         .169395306766867743169300202489,
                         .380690406958401545684749139159,
                         .619309593041598454315250860841,
                         .830604693233132256830699797511,
                         .966234757101576013906150777244
        };
        DataType w6[] = {
                         .085662246189585172520148071111,
                         .180380786524069303784916756904,
                         .233956967286345523694935171995,
                         .233956967286345523694935171995,
                         .180380786524069303784916756904,
                         .085662246189585172520148071111
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x6, x6 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w6, w6 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 7;

        DataType x7[] = {
                         .025446043828620737736905157977,
                         .12923440720030278006806761336,
                         .297077424311301416546696793961,
                         .5,
                         .702922575688698583453303206039,
                         .87076559279969721993193238664,
                         .974553956171379262263094842023
        };
        DataType w7[] = {
                         .064742483084434846635305716358,
                         .139852695744638333950733885707,
                         .190915025252559472475184887747,
                         .208979591836734693877551020408,
                         .190915025252559472475184887747,
                         .139852695744638333950733885707,
                         .064742483084434846635305716358
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x7, x7 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w7, w7 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 8;

        DataType x8[] = {
                         .019855071751231884158219565716,
                         .101666761293186630204223031762,
                         .237233795041835507091130475406,
                         .408282678752175097530261928822,
                         .591717321247824902469738071178,
                         .762766204958164492908869524594,
                         .898333238706813369795776968238,
                         .980144928248768115841780434284
        };
        DataType w8[] = {
                         .050614268145188129576265677161,
                         .111190517226687235272177997216,
                         .156853322938943643668981100988,
                         .18134189168918099148257522464,
                         .18134189168918099148257522464,
                         .156853322938943643668981100988,
                         .111190517226687235272177997216,
                         .050614268145188129576265677161
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x8, x8 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w8, w8 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 9;

        DataType x9[] = {
                         .015919880246186955082211898548,
                         .081984446336682102850285105965,
                         .193314283649704801345648980329,
                         .337873288298095535480730992679,
                         .5,
                         .662126711701904464519269007321,
                         .806685716350295198654351019671,
                         .918015553663317897149714894035,
                         .984080119753813044917788101452
        };
        DataType w9[] = {
                         .040637194180787205985946079062,
                         .090324080347428702029236015625,
                         .130305348201467731159371434718,
                         .156173538520001420034315203291,
                         .165119677500629881582262534644,
                         .156173538520001420034315203291,
                         .130305348201467731159371434718,
                         .090324080347428702029236015625,
                         .040637194180787205985946079062
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x9, x9 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w9, w9 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 10;

        DataType x10[] = {
                          .013046735741414139961017993958,
                          .067468316655507744633951655788,
                          .160295215850487796882836317442,
                          .283302302935376404600367028417,
                          .425562830509184394557586999435,
                          .574437169490815605442413000565,
                          .716697697064623595399632971583,
                          .839704784149512203117163682558,
                          .932531683344492255366048344212,
                          .986953264258585860038982006042
        };
        DataType w10[] = {
                          .0333356721543440687967844048950,
                          .0747256745752902965728881698850,
                          .109543181257991021997767467126,
                          .134633359654998177545613460784,
                          .147762112357376435086946497326,
                          .147762112357376435086946497326,
                          .134633359654998177545613460784,
                          .109543181257991021997767467126,
                          .0747256745752902965728881698850,
                          .0333356721543440687967844048950
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x10, x10 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w10, w10 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 11;

        DataType x11[] = {
                          .010885670926971503598030999438,
                          .056468700115952350462421115348,
                          .134923997212975337953291873984,
                          .240451935396594092037137165270,
                          .365228422023827513834234007299,
                          .5,
                          .634771577976172486165765992701,
                          .759548064603405907962862834730,
                          .865076002787024662046708126016,
                          .943531299884047649537578884652,
                          .989114329073028496401969000562
        };
        DataType w11[] = {
                          .0278342835580868332413768601148,
                          .0627901847324523123173471498315,
                          .0931451054638671257130488207810,
                          .116596882295995239959261852412,
                          .131402272255123331090344434948,
                          .136462543388950315357241764168,
                          .131402272255123331090344434948,
                          .116596882295995239959261852412,
                          .0931451054638671257130488207810,
                          .0627901847324523123173471498315,
                          .0278342835580868332413768601148
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x11, x11 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w11, w11 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 12;

        DataType x12[] = {
                          .009219682876640374654725454925,
                          .047941371814762571660767066940,
                          .115048662902847656481553083393,
                          .206341022856691276351648790529,
                          .316084250500909903123654231678,
                          .437383295744265542263779315268,
                          .562616704255734457736220684732,
                          .683915749499090096876345768322,
                          .793658977143308723648351209471,
                          .884951337097152343518446916607,
                          .952058628185237428339232933060,
                          .990780317123359625345274545075
        };
        DataType w12[] = {
                          .0235876681932559135973079808979,
                          .0534696629976592154801273590260,
                          .0800391642716731131673262646270,
                          .101583713361532960874532227889,
                          .116746268269177404380424949464,
                          .124573522906701392500281218021,
                          .124573522906701392500281218021,
                          .116746268269177404380424949464,
                          .101583713361532960874532227889,
                          .0800391642716731131673262646270,
                          .0534696629976592154801273590260,
                          .0235876681932559135973079808979
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x12, x12 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w12, w12 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 13;

        DataType x13[] = {
                          .007908472640705925263585275596,
                          .041200800388511017396726081749,
                          .099210954633345043602896755208,
                          .178825330279829889678007696502,
                          .275753624481776573561043573936,
                          .384770842022432602967235939451,
                          .5,
                          .615229157977567397032764060549,
                          .724246375518223426438956426064,
                          .821174669720170110321992303498,
                          .900789045366654956397103244792,
                          .958799199611488982603273918251,
                          .992091527359294074736414724404
        };
        DataType w13[] = {
                          .0202420023826579397600107964706,
                          .0460607499188642239572108885002,
                          .0694367551098936192318008890105,
                          .0890729903809728691400233461840,
                          .103908023768444251156261609634,
                          .113141590131448619206045093020,
                          .116275776615436955097294757635,
                          .113141590131448619206045093020,
                          .103908023768444251156261609634,
                          .0890729903809728691400233461840,
                          .0694367551098936192318008890105,
                          .0460607499188642239572108885002,
                          .0202420023826579397600107964706
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x13, x13 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w13, w13 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 14;

        DataType x14[] = {
                          .006858095651593830579201366648,
                          .035782558168213241331804430311,
                          .086399342465117503405102628675,
                          .156353547594157264925990098490,
                          .242375681820922954017354640724,
                          .340443815536055119782164087916,
                          .445972525646328168966877674890,
                          .554027474353671831033122325110,
                          .659556184463944880217835912084,
                          .757624318179077045982645359276,
                          .843646452405842735074009901510,
                          .913600657534882496594897371325,
                          .964217441831786758668195569689,
                          .993141904348406169420798633352
        };
        DataType w14[] = {
                          .0175597301658759315159164387279,
                          .0400790435798801049028166394764,
                          .0607592853439515923447074066035,
                          .0786015835790967672848009696255,
                          .0927691987389689068708582949900,
                          .102599231860647801982962032827,
                          .107631926731578895097938221659,
                          .107631926731578895097938221659,
                          .102599231860647801982962032827,
                          .0927691987389689068708582949900,
                          .0786015835790967672848009696255,
                          .0607592853439515923447074066035,
                          .0400790435798801049028166394764,
                          .0175597301658759315159164387279
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x14, x14 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w14, w14 + size );

        index_field_[size] = cntr;
        cntr++;

        // ---------------------------------------------
        size = 15;

        DataType x15[] = {
                          .006003740989757285755217140706,
                          .031363303799647047846120526145,
                          .075896708294786391899675839613,
                          .137791134319914976291906972693,
                          .214513913695730576231386631373,
                          .302924326461218315051396314509,
                          .399402953001282738849685848302,
                          .5,
                          .600597046998717261150314151698,
                          .697075673538781684948603685491,
                          .785486086304269423768613368627,
                          .862208865680085023708093027307,
                          .924103291705213608100324160387,
                          .968636696200352952153879473855,
                          .993996259010242714244782859294
        };
        DataType w15[] = {
                          .0153766209980586341773141960552,
                          .0351830237440540623546337059413,
                          .0535796102335859675059347752550,
                          .0697853389630771572239023974620,
                          .0831346029084969667766004305505,
                          .0930805000077811055134002809275,
                          .0992157426635557882280591632215,
                          .101289120962780636440310099984,
                          .0992157426635557882280591632215,
                          .0930805000077811055134002809275,
                          .0831346029084969667766004305505,
                          .0697853389630771572239023974620,
                          .0535796102335859675059347752550,
                          .0351830237440540623546337059413,
                          .0153766209980586341773141960552
        };

        x_[cntr].reserve ( size );
        x_[cntr].insert ( x_[cntr].begin ( ), x15, x15 + size );
        y_[cntr].resize ( size, 0.0 );
        z_[cntr].resize ( size, 0.0 );
        w_[cntr].reserve ( size );
        w_[cntr].insert ( w_[cntr].begin ( ), w15, w15 + size );

        index_field_[size] = cntr;
        cntr++;

    }

    template<class DataType>
    QuadratureGaussLine<DataType>* QuadratureGaussLine<DataType>::clone ( ) const
    {
        return new QuadratureGaussLine<DataType>( *this );
    }

    // template instanciation
    template class QuadratureGaussLine<double>;
    template class QuadratureGaussLine<float>;

} // namespace hiflow
