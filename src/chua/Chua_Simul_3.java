
package chua;

/*
 * this is file : \Documents\NetBeansProjects\ChuaOscillator\src\chua\Chua_Simul_3.java.
 * simulate a Poincare map of degree 3 in the x'-y' plane
 * based on standard cubic model, Book IV, page 62 and page 69.

 * separately: convert real, 2D, cubic, coeff Cxy to complex gij - July, 2024
 * multiply complex gij to produce 2D point
 * implement, and invert, a cubic transform in z-space, using gij
 * to emulate a normalized cubic model of a torus
 * 'hdr' must be the output of 'Chua_N_response_circ.py'
 * Note: for Cxy we use a series expansion in (x, y) without any factorials
 * Note: for gij, hij, we use standard Taylor expansion with factorials
 */

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.*;
import javax.swing.*;

public class Chua_Simul_3 extends JDialog
{
    //private static String hdr = " 69597600 ,  99.9948 ,  1499.25037 ,  -0.51325 ,  -1.0 ,  0.144 ,  2400 ,  6.990240639707452E-5 , 1.2809579565997563e-06 , 0.9826226599215292 , -0.18659758937064014 , -1.103029973043194 , 8.462874528763841 , -34.24146413987307 , 8.728918049582369 , 272.23195311295456 , -1730.6613177460904 , -1879.1440596659322 , -1.0364954821513763e-06 , 0.18683689988920496 , 0.9821346041621322 , -0.004294046470095772 , 43.389922533105086 , 9.147510786513596 , 44.33005570148138 , 1047.6157267475264 , -440.64698522425294 , 3090.601610076812";
    //private static String hdr = " 69597600 ,  99.9948 ,  1499.25037 ,  -0.51325 ,  -1.0 ,  0.144 ,  2400 ,  6.990240639707452E-5 , 0 , 0.9825702656507906 , -0.18687170000740827 , -1.0745125402430635 , 8.996841269738953 , -32.42897400057778 , 9.369569463824062 , 232.6286221161031 , -1698.9717975640954 , -1868.6516702490385 , 0 , 0 , 0 , 0 , 0 , 0 , 0.18687170000740827 , 0.9825702656507906 , 0.08759955917253591 , 42.46383864701076 , 7.4889898795593615 , 44.529932224129105 , 1034.154272302978 , -415.0566709553265 , 3083.976801672329 , 0 , 0 , 0 , 0 , 0";
    //private static String hdr = " 347997600 ,  99.98 ,  1499.25037 ,  -0.51325 ,  -1.0 ,  0.144 ,  2400 ,  6.990205973471885E-5 , 0 , 0.9825013928586689 , -0.18646770776055427 , -1.2306352795346436 , 9.45266581283182 , -32.8516447072305 , 4.932463161858538 , 278.00455317334763 , -1642.4571438790335 , -2005.6936708215628 , -77.2460484547081 , 4871.269477570737 , -60310.748588801536 , -63415.74647710711 , -226428.16048929538 , 0 , 0.18646770776055427 , 0.9825013928586689 , -0.6151693710061411 , 42.91447086704597 , 7.701425537642308 , 28.50130065842322 , 1083.4982863191667 , -576.6878374059502 , 3137.8062713894487 , 1207.5983017416697 , 17839.211077886783 , -59774.20457562876 , 227869.81960854883 , 216790.45046234503";
    //private static String hdr = " 0 ,  0.6155 ,  0.6 ,  1.25 ,  NaN ,  NaN ,  2400 ,  -0.00203955745286074 , 0 , 0.9677529263313476 , -0.2574846644096744 , 1.2147694074060806 , -4.76804647344731 , 6.487280766484504 , 12.286191019464926 , -66.1426052057399 , 113.83491403046077 , -81.79561452649884 , -1855.550090916026 , 383.161814218651 , 685.6070150780075 , 534.9640154506382 , 1460.214227204807 , 0 , 0.2574846644096744 , 0.9677529263313476 , -1.522974464099001 , 2.1102407238384404 , -1.2024838653488197 , -2.2001957485654806 , -2.372802246687878 , 11.884344420089013 , -11.425791166636012 , -1729.0849338313662 , 798.1907030734329 , -457.4790199136309 , 1685.4357633940065 , 685.0745241195965";
    //private static String hdr = " 347997700 ,  99.98 ,  1499.25037 ,  -0.51325 ,  -1.0 ,  0.144 ,  2400 ,  6.990205973471885E-5 , 0 , 0.9825013928233131 , -0.18646770775949692 , -3.98 , 38.3 , -31.5 , 2650.0 , -6684.0 , 6738.0 , -2616.0 , 0 , 0 , 0 , 0 , 0 , 0 , 0.18646770775949692 , 0.9825013928233131 , -36.35 , 18.13 , 9.73 , -1860 , 870 , 959 , -176 , 0 , 0 , 0 , 0 , 0";
    //private static String hdr = " 997920000 ,  99.986 ,  1499.25037 ,  -0.51325 ,  -1.0 ,  0.144 ,  2400 ,  6.99022002895545E-5 , 0 , 0.9825294302676153 , -0.18663142507958805 , -1.258844279421771 , 9.208782516253333 , -34.726606220615984 , 4.926786225224641 , 277.9084138764038 , -1642.598947865595 , -2008.3910011695327 ,0,0,0,0,0, 0 , 0.18663142507958805 , 0.9825294302676153 , -0.7302370019166385 , 43.91503096442345 , 9.357872279018837 , 28.487868661445887 , 1083.7092119043655 , -575.655997201535 , 3139.7502773965234,0,0,0,0,0";
    //private static String hdr = " 1396320000 ,  99.983 ,  1499.25037 ,  -0.51325 ,  -1.0 ,  0.144 ,  2400 ,  6.990212988270638E-5 , 0 , 0.9825146336901965 , -0.18654952002889488 , -1.233689396260997 , 9.488325814359307 , -32.68445419965168 , 3.73856987429206 , 262.0202237095138 , -1531.3319874178137 , -1786.4871886749431 , 205.18854428134958 , 5373.6630429863035 , -62559.901949079685 , -74235.70484604409 , -249350.8029182552 , 25365.55698296425 , 36678.33005782716 , -1724157.2772215619 , -333549.5908996819 , -16904891.503100708 , -22130456.078911066 , 0 , 0.18654952002889488 , 0.9825146336901965 , -0.6071492854165632 , 42.838486758943006 , 7.464786773514292 , 33.17501121603761 , 1073.4816667823288 , -661.7061541277534 , 2805.5815922987713 , 960.0895382158386 , 18339.620871800427 , -60857.409575948426 , 248004.51393209893 , 250358.18237968118 , 9704.459382312793 , 44051.765264390306 , -2745510.953070757 , 11167658.619492525 , 15024219.720379766 , 32505996.41001745";

    //private static String hdr = " 674880000 ,  99.994 ,  1499.25037 ,  -0.51325 ,  -1.0 ,  0.144 ,  240000 ,  6.990238766388633E-7 , 0 , 0.9825665599529959 , -0.1868499075617937 , -5.813871391585502 , 41.7895205462587 , -12.03812472182267 , 3322.0127631249993 , -2172.661960459033 , 1348.0255431749335 , -334.69805591279743 , -168843.58733845723 , 445861.49092453445 , -248132.66901944793 , 63623.625743015116 , -6989.328603287862 , 24383896.867314197 , -43459424.2592648 , 41988310.67266308 , -19984072.880590297 , 3594459.845558569 , -410408.87236746424 , 0 , 0.1868499075617937 , 0.9825665599529959 , -34.1222051411813 , -0.9366988711967978 , 2.2094969605497417 , 480.00423631029895 , -2662.6481336618926 , 862.1680986529439 , -16.44550486346501 , -165105.54721105774 , 156750.09594828307 , -143079.56754445593 , 56024.31695680063 , -4393.189171842834 , 8388270.466699365 , -32076918.546389773 , 8133093.290679754 , -21985422.66464864 , -5716678.054900305 , -2382468.5945130163";
    //private static String hdr = " 0, 0.6137, 0.6, 1.25, NaN, NaN, NaN, -0.0010 , 0 , 0.9815763708664 , -0.2756924264915933, -0.04147536672454766, 0.49992905987251607, 1.628094853959384, 0.34778776015627666, -0.3757385191356548, -1.8345052644144124, -1.129625223257578 , 0 , 0.2756924264915933 , 0.9815763708664, -0.6851662571560203, -1.1360709075611275, -2.698598016167229, 1.2562576979100823, 0.9573019393852975, 2.4385172774278927, 11.599658895680303";
    //private static String hdr = " 73920000, 99.986, 1499.25037, -0.51325, -1.0, 0.144, NaN, 5.0E-5 , 0 , 0.982529316888481 , -0.18663155943604884, 1.4979171715092254, 41.42974807788249, -6.417687588934195, 3394.3715374118374, -3543.6824421438187, 1159.126459445238, -1280.4601332355646, -709052.0076599567, 766346.1684626958, -1029818.4769828893, 68564.88436072569, -280488.57464032713 , 0 , 0.18663155943604884 , 0.982529316888481, -28.75291067582463, -1.261018152589866, 5.638430326813317, 340.06133400046355, -3601.2606807769944, 775.9330414509485, -596.6475894573047, -537670.9378060552, 370728.66898158984, -673718.3820685825, 63634.93393505788, -171674.35925809544";
    //private static String hdr = " 0, 0.615, 0.6, 1.25, NaN, NaN, NaN, -0.0010 , 0 , 0.9715853168026777 , -0.26265003717794533, 0.033836616009405865, 0.5341282853889122, 1.74555962734054, -0.46624591949623795, -0.08252239434874711, -2.9847430716242758, -2.6083242241060467, 11.686416376791614, -11.826295576551024, 27.383780056585593, -6.0806677417808865, 11.638505282228962 , 0 , 0.26265003717794533 , 0.9715853168026777, -0.5647038441605761, -1.1764071266333436, -2.4560462162720276, -0.007757200854699303, 3.5125081898908235, 3.6034816708189696, 10.979341843952149, 20.10323909857962, -25.746116566238957, 58.82991591440729, -56.95374513727638, 6.069732680592089";
    //private static String hdr = " 73920000 , 99.986 , 1499.25037 , -0.51325 , -1.0 , 0.144 , NaN , 5.0E-5  , 0, 0.982529316888481 , -0.18663155943604884 , 0, 0, 0, 3450.290775309426 , -13437.034867247576 , 1064.520478988718 , 6026.099432700741 , -940751.8753374612 , -242901.12527098507 , 1285921.1587427037 , 2241615.589091816 , 320865.4461114729 ,40601344.878941864 , 456465865.7944363 , 9652458613.27098 , 3540949643.345993 , -12635625332.298664 , -1895284355.563605 , 0, 0.18663155943604884 , 0.982529316888481 , 0, 0, 0, 2807.5634338954133 , -1448.0820946060435 , -10687.764356793112 , -1278.498360652526 , 17303.272907657527 , 1378705.3642015634 , -4682151.813397732 , -832721.8086269677 , 370159.32608726626 ,26438309.22239845 , -3167515265.4028573 , -710832704.8221453 , 16297983518.285196 , 4272479189.940985 , -2949341006.523123";
    //private static String hdr = " 73920000, 99.986, 1499.25037, -0.51325, -1.0, 0.144, NaN, 5.0E-5 , 0 , 0.982529316888481 , -0.18663155943604884, 1.4979171715092254, 41.42974807788249, -6.417687588934195, 3394.3715374118374, -3543.6824421438187, 1159.126459445238, -1280.4601332355646, -709052.0076599567, 766346.1684626958, -1029818.4769828893, 68564.88436072569, -280488.57464032713 , 0 , 0.18663155943604884 , 0.982529316888481, -28.75291067582463, -1.261018152589866, 5.638430326813317, 340.06133400046355, -3601.2606807769944, 775.9330414509485, -596.6475894573047, -537670.9378060552, 370728.66898158984, -673718.3820685825, 63634.93393505788, -171674.35925809544";

    //private static String hdr = " 69597600  ,  99.9948  ,  1499.25037  ,  -0.51325  ,  -1.0  ,  0.144  ,  2400  ,  6.990240639707452E-5  , 0, 0.9825702656507906 , -0.18687170000740827 , 0, 0, 0, -1457.1240443472166 , 9782.854979882632 , 7741.359223003977 , -242.9945038075375 , 202609.08531924285 , 1935276.8283106992 , -2317595.0568399434 , -1703618.4712000187 , 38595.0859715058 ,95162173.77159268 , -182055979.18801197 , -1216906389.105405 , 219587980.0958921 , 338043683.1029007 , -7195850.349761835 , 0, 0.18687170000740827 , 0.9825702656507906 , 0, 0, 0, -7672.315290850676 , -10019.5307173032 , 5922.370247765803 , 944.1795584454187 , -413378.9317032323 , 2212957.822224167 , 3197956.396671033 , -21407.65200155339 , -29691.287421174842 ,37464400.47728396 , 942127581.3069272 , -585121659.0199361 , -796992655.5638853 , -33535773.044681508 , 1950994.8251135598";
    //private static String hdr = " 564960000, 99.992, 1499.25037, -0.51325, -1.0, 0.144, NaN, 5.0E-5 , 0 , 0.9825572305065098 , -0.18679531980262137, -4.815319687213514, 39.57022741061535, -9.826177723762106, 4192.354649346585, -2937.0258675850237, 1344.2018973090053, -650.6168522102133, -372010.38343923836, 753680.6526983386, -422373.1423822457, 71457.68178773094, -51048.05810241968 , 0 , 0.18679531980262137 , 0.9825572305065098, -32.28579234404382, -1.807796595288989, 3.3896307248343227, 653.0360027361866, -3302.3011498791643, 885.5337739120685, -187.59514996935627, -299131.1124235608, 335234.0452514214, -279791.57281488826, 54674.50161004579, -27585.285916760957";

    //private static String hdr = " 73920000 , 99.986 , 1499.25037 , -0.51325 , -1.0 , 0.144 , NaN , 5.0E-5  , 0, 0.982529316888481 , -0.18663155943604884 , 0, 0, 0, 3450.277338133448 , -13437.082332153232 , 1065.0826780286543 , 6025.685344572169 , -940844.3341580641 , -242781.1826756078 , 1285857.5837793644 , 2241523.297230466 , 320899.9780091982 ,60085589.606249005 , -288342958.8691498 , -414999582.25601375 , 2309053782.789438 , 418261529.36426926 , -174918245.73308337 , 0, 0.18663155943604884 , 0.982529316888481 , 0, 0, 0, 2807.5453077913726 , -1448.2871249919463 , -10687.373069690195 , -1278.5713310937513 , 17195.79193606115 , 1378805.305077914 , -4682143.9084503995 , -832550.3358831485 , 370120.8373641219 ,75477750.74504101 , -24548107.856253743 , -1029019632.426597 , 366730101.21108425 , 413448447.0426211 , 8014223.735627143";
    //private static String hdr = " 748317600, 99.98, 1499.25037, -0.51325, -1.0, 0.144, NaN, 5.0E-5 , 0 , 0.9825013613218427 , -0.18646772234613535, 4.042963542757758, 0.810381882305268, -7.805917846886524, 63.27972656500435, 125.09584444463141, 85.93186315497067, 69.89125424586992, 597.1733040620396, 542.2933886275923, -277.68284126536344, -2116.4807448993847, -728.7754358053828 , 0 , 0.18646772234613535 , 0.9825013613218427, 6.258785764051628, 6.73889115634888, -1.4952981764436284, 86.58429528379921, 126.93930384081399, 5.4481185881293115, 34.20118805297246, 601.6321889245895, 1317.788059517367, -523.8084786233038, -454.1043331199748, -182.15097917393376";
    //private static String hdr = " 748319900 , 99.98 , 1499.25037 , -0.51325 , -1.0 , 0.144 , NaN , 5.0E-5  , 0, 0.9825013637459529 , -0.1864677328353816 , 0, 0, 0, -274.84333944383843 , 1626.256033504982 , 435.2552194705163 , -261.9583681881735 , 70447.69139543881 , -143438.94124141027 , -41147.32886741769 , 13014.647816699493 , 1008.4859079515372 ,4021469.9669092814 , -7904419.876289886 , 4689834.77870477 , 21540489.558443524 , -10967481.91138449 , 1055756.831287527 , 0, 0.1864677328353816 , 0.9825013637459529 , 0, 0, 0, -1178.2240966922666 , 2302.4428266351642 , 780.9747104094367 , -220.56159737221788 , -41885.395972904815 , 98423.67180481146 , 42341.56169759987 , -74823.47139161959 , 18769.991146315133 ,11626389.207679009 , -30611002.16687298 , 3916319.834857713 , 2573423.1440527067 , 2193661.5839719023 , -972495.4797072884";
    //private static String hdr = " 1026720000, 99.978, 1499.25037, -0.51325, -1.0, 0.144, NaN, 5.0E-5 , 0 , 0.9824922357197969 , -0.18641319289687444, -16.794633630773355, 83.63543404276152, -14.308159197240348, -17332.38966393748, 2359.601928388419, -10612.990034704204, 15453.064869495916, -4857060.857610761, -4799743.6571988445, -3963818.6139985393, -5972961.924411841, 1109704.8579936635 , 0 , 0.18641319289687444 , 0.9824922357197969, 70.89322226858978, 20.01857282716222, 24.870067962083606, 81452.58352873953, 15040.936618239051, 82441.77221822333, 23375.798677093793, 20303871.437058974, -81242.92969483383, 41585903.45680854, 1797869.8237645745, 21001970.459979594";
    //private static String hdr = " 115200000, 99.9948, 1499.25037, -0.51325, -1.0, 0.144, NaN, 5.0E-5 , 0 , 0.9825695507983214 , -0.1868718579752067, 17.1152105256668, 46.36643415540755, -38.59138476769669, 4.84998395177672, -1863.9219567978432, 2580.6365868428807, 373.29179463545177, -27838.08595845153, -294503.5264571201, 220.99607492865744, 141143.6095737989, 32662.31164886231 , 0 , 0.1868718579752067 , 0.9825695507983214, -64.58576956648034, 7.734256106161378, 40.15093074371333, -3710.4799255599933, 6476.2762330665255, -8892.81553580754, 4671.314640950858, -58324.72689410096, 409633.4782953956, -388489.408643748, -17099.06981794391, 83259.07052468101";
    //private static String hdr = "0 , 2.1 , NaN, NaN, NaN, NaN, NaN, NaN , 0 , 0.5 , -0.9219544457292888, 0, 0, 0, 2.804184010263312, -3.096053990269489, 1.7280138593297225, 2.3457639618867354 , 0 , 0.9219544457292888 , 0.5, 0, 0, 0, -2.500494394199581, 3.2772226907500546, -1.4243242432659944, -2.1645952614061734";

    //private static String hdr = "13 , 2.20 , NaN, NaN, NaN, NaN, NaN, NaN , 0 , 0.49999999999999994 , -0.9746794344808964, 0, 0, 0, 2.5031494545806052, -3.1401971979326717, 1.7474112105536945, 2.295439264923758 , 0 , 0.9746794344808964 , 0.49999999999999994, 0, 0, 0, -2.286230665914016, 3.2737295867321428, -1.5304924218871094, -2.161906876124286";
    //private static String hdr = "0 , 2.14 , NaN, NaN, NaN, NaN, NaN, NaN , 0 , 0.49999999999999994 , -0.9433981132056605, -1.6370988658034027, -1.134197731606805, 0.5029011341965977, 0, 0, 0, 0 , 0 , 0.9433981132056605 , 0.49999999999999994, 1.6370988658034025, 1.1341977316068048, -0.5029011341965975, 0, 0, 0, 0";

    //private static String hdr = "0 , -0.36 , 1.024 , -0.1 , NaN, NaN, NaN, NaN , 0 , 0.20946988283227083 , -0.9798841031140071, -0.782489631991973, 0.9137287079086109, -0.2635498663274307, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0.9798841031140071 , 0.20946988283227083, 0.782489631991973, -0.9137287079086109, 0.2635498663274307, 0, 0, 0, 0, 0, 0, 0, 0, 0";
    //private static String hdr = "0 , -0.36 , 1.0249777 , -0.1 , NaN, NaN, NaN, NaN , 0 , 0.209340771971716 , -0.9804167097780551, -0.7825746941328198, 0.9143655752339029, -0.2638931501893879, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0.9804167097780551 , 0.209340771971716, 0.7825746941328199, -0.914365575233903, 0.26389315018938797, 0, 0, 0, 0, 0, 0, 0, 0, 0";
    //private static String hdr = "0 , -0.2 , 1.018 , -0.1 , NaN, NaN, NaN, NaN , 0 , 0.11038956094277855 , -0.9976476599062921, -0.670853177261784, 0.9743680871745739, -0.350074055376227, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0.9976476599062921 , 0.11038956094277855, 0.670853177261784, -0.9743680871745739, 0.350074055376227, 0, 0, 0, 0, 0, 0, 0, 0, 0";
    //private static String hdr = "0 , 2.22 , NaN, NaN, NaN, NaN, NaN, NaN , 0 , 0.5 , -0.9848857801796106, -1.6735171216491587, -1.127034243298317, 0.5464828783508415, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0.9848857801796106 , 0.5, 1.6735171216491587, 1.127034243298317, -0.5464828783508415, 0, 0, 0, 0, 0, 0, 0, 0, 0";
    //private static String hdr = "0 , 2.14 , NaN, NaN, NaN, NaN, NaN, NaN , 0 , 0.49999999999999994 , -0.9433981132056605, -1.6370988658034027, -1.134197731606805, 0.5029011341965977, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0.9433981132056605 , 0.49999999999999994, 1.6370988658034025, 1.1341977316068048, -0.5029011341965975, 0, 0, 0, 0, 0, 0, 0, 0, 0";
    //private static String hdr = "0 , 2.1 , NaN, NaN, NaN, NaN, NaN, NaN , 0 , 0.5 , -0.9219544457292888, -1.6194424517739725, -1.138884903547945, 0.4805575482260276, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0.9219544457292888 , 0.5, 1.6194424517739725, 1.138884903547945, -0.4805575482260276, 0, 0, 0, 0, 0, 0, 0, 0, 0";

    //private static String hdr = "0   , 2.08   , NaN  , NaN  , NaN  , NaN  , NaN  , NaN   , 0, 0.49999999999999994 , -0.9110433579144299 , 0, 0, 0, 0.8901528484684819 , 1.6301434605614837 , 0.8901528484684819 , 1.630143460561484 , 0 , 0.0 , 0 , 0.0 , 0 ,0 , 0.0 , 0 , 0.0 , 0 , 0.0 , 0, 0.9110433579144299 , 0.49999999999999994 , 0, 0, 0, -1.6301434605614844 , 0.8901528484684814 , -1.6301434605614846 , 0.8901528484684814 , 0 , 0.0 , 0 , 0.0 , 0 ,0 , 0.0 , 0 , 0.0 , 0 , 0.0";
    //private static String hdr = "2  , -0.36  , 1.022  , -0.1  , NaN , NaN , NaN , NaN  , 0, 0.20973455993273862 , -0.9787934891366136 , 0, 0, 0, 0.4613732764545197 , 0.11673926950595609 , 0.4613732764545197 , 0.11673926950595609 , 0 , 0.0 , 0 , 0.0 , 0 ,0 , 0.0 , 0 , 0.0 , 0 , 0.0 , 0, 0.9787934891366136 , 0.20973455993273862 , 0, 0, 0, -0.11673926950595609 , 0.4613732764545197 , -0.11673926950595609 , 0.4613732764545197 , 0 , 0.0 , 0 , 0.0 , 0 ,0 , 0.0 , 0 , 0.0 , 0 , 0.0";
// Dec 8, 2025 test cubic theory
    //private static String hdr = "0 , 2.1 , NaN, NaN, NaN, NaN, NaN, NaN , 0 , 0.5 , -0.9219544457292888, 0, 0, 0, 2.742923552471728, -3.338483148541635, 1.666753401538138, 2.1033348036145885 , 0 , 0.9219544457292888 , 0.5, 0, 0, 0, -2.7429235524717273, 3.3384831485416386, -1.6667534015381404, -2.1033348036145885";
    //private static String hdr = "0 , 2.223 , NaN, NaN, NaN, NaN, NaN, NaN , 0 , 0.5 , -0.99498743710662, -1.6828211765451586, -1.1256423530903172, 0.5571788234548416, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0.99498743710662 , 0.5, 1.6828211765451586, 1.1256423530903172, -0.5571788234548416, 0, 0, 0, 0, 0, 0, 0, 0, 0";
    //private static String hdr = "0 , 2.174 , NaN, NaN, NaN, NaN, NaN, NaN , 0 , 0.5 , -0.9612491872558333, -1.6524101009453953, -1.1308202018907907, 0.5215898990546046, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0.9612491872558333 , 0.5, 1.6524101009453958, 1.130820201890791, -0.5215898990546047, 0, 0, 0, 0, 0, 0, 0, 0, 0";
    private static String hdr = "0 , 2.27 , NaN, NaN, NaN, NaN, NaN, NaN , 0 , 0.5 , -0.9628083921528728, -1.653773427647466, -1.1305468552949316, 0.5232265723525342, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 0 , 0.9628083921528728 , 0.5, 1.653773427647466, 1.1305468552949316, -0.5232265723525342, 0, 0, 0, 0, 0, 0, 0, 0, 0";

    private static double[][] C_g_2_R = new double[][] {{ 1, 0,-1},     // convert C to g
                                                        { 1, 0, 1},     // quadratic, real
                                                        { 1, 0,-1}};
    private static double[][] C_g_2_I = new double[][] {{ 0,-1, 0},     // see 'Chua_2D_cubic_variable_c1.py'
                                                        { 0, 0, 0},
                                                        { 0, 1, 0}};
    private static double[][] C_g_3_R = new double[][] {{ 3, 0,-3, 0},
                                                        { 3, 0, 1, 0},
                                                        { 3, 0, 1, 0},
                                                        { 3, 0,-3, 0}};
    private static double[][] C_g_3_I = new double[][] {{ 0,-3, 0, 3},
                                                        { 0,-1, 0,-3},
                                                        { 0, 1, 0, 3},
                                                        { 0, 3, 0,-3}};
    private static Point2D.Double[] h_2 = new Point2D.Double[3];        // quadratic response
    private static Point2D.Double[] inv_3 = new Point2D.Double[4];      // inverse cubic response

    private static final BufferedImage image = new BufferedImage(480, 480, BufferedImage.TYPE_3BYTE_BGR);
    private static final Graphics2D DC = image.createGraphics();
    private static final JLabel lblImage = new JLabel(new ImageIcon(image));
    private static JButton btnCalc = new JButton("Calc");
    private static JButton btnResetHenon = new JButton("Res_Henon");
    private static JButton btnReset = new JButton("Reset");
    private static JButton btnClear = new JButton("Clear");
    private static JButton btnInc = new JButton("+");
    private static JCheckBox invertChk = new JCheckBox(" invert ");
    private static JCheckBox printChk = new JCheckBox("  print  ");

    private static JTextField txtalpha = new JTextField();
    private static JTextField txttheta = new JTextField();
    private static JTextField txta = new JTextField();
    private static JTextField txtb = new JTextField();
    private static JTextField[][] txtCarr = {{new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()},
                                             {new JTextField(), new JTextField()}};
    private static String Chua_iTstr;
    private static double Chua_alpha, Chua_beta, Chua_gamma, Chua_a, Chua_c, Chua_Period, Chua_delt;
    private static double[][] Carr;
    private static JTextField txtstartx = new JTextField("0.0"); //  0.001 // -0.2
    private static JTextField txtstarty = new JTextField("0.0"); // -0.001 // -0.2
    private static JTextField txtrange = new JTextField("1");   // 0.005
    private static JTextField txtbeta_override = new JTextField("0");
    private static JTextField txtcolor_override = new JTextField("0");
    private static double x0, y0;
    //private static double xclick, yclick;
    private static double Henon_a = -0.36; // -0.050; // 3; // 1.89; // -0.20; //-0.36; //-0.21; // -0.2285; // 0.3;
    private static double Henon_b = 1.023; // 1.0249777379; //1.024977738; // 1.0146; // 1.495; // 0.9323; // 1.02496; //-1.055; //-1.057; // 1.02; //1.011; //1.023; // 1.024977; // 1.0249778;  // 1.0106;
    private static double Henon_R = -0.1; // 0.5; //-0.50; // -0.1; // -0.5;
    private static double Henon_S = 0;
    //private static double Delay_a = 2.2553; // 2.2559; // 2.272; // 2.1; // 2.2027; // 2.230; // 2.271;
    //private static final double H_cos = 0.00; // 0.00; // -0.45; // 0.24; // 0.8;           // Henon map from Arrowsmith p.54
    //private static double H_alpha = Math.acos(H_cos);   // Henon map from Arrowsmith p.54
    private static int iT = 0;

    public Chua_Simul_3()
    {
        //setTitle("Chua System - Simulate x'-y' Scatter (quartic map) (" + hdr.split(",")[0].trim() + "_" + hdr.split(",")[1].trim() + "_" + (hdr.split(",").length/2 - 7) + ")");
        //setTitle("Delayed Logistic - Simulate x-y Scatter : a = " + Delay_a);
        setTitle("Henon_scatter: (" + Henon_a + ", " + Henon_b + ", " + Henon_R + ", " + Henon_S + ")");
        //setTitle("Henon_1969_scatter: (" + H_cos + ", " + String.format("%.8f", 2*Math.tan(H_alpha/2)) + ", " + String.format("%.8f", 2*Math.tan(H_alpha/2)*Math.tan(H_alpha/2)) + ")");
        setIconImage(Toolkit.getDefaultToolkit().getImage(Main.class.getResource("images/icon.gif")));
        setSize(780, 544);
        setLocationByPlatform(true);

        DC.setBackground(Color.white);
        DC.clearRect(0, 0, image.getWidth(), image.getHeight());
        lblImage.setBorder(BorderFactory.createEtchedBorder());

        final JLabel[] lblarr = {new JLabel("C 20"),
                                 new JLabel("C 11"),
                                 new JLabel("C 02"),
                                 new JLabel("C 30"),
                                 new JLabel("C 21"),
                                 new JLabel("C 12"),
                                 new JLabel("C 03")};
        JPanel[] spacerPanel = new JPanel[2];
        for (int i = 0; i < spacerPanel.length; i++)
        {
            spacerPanel[i] = new JPanel();
            spacerPanel[i].setPreferredSize(new Dimension(190, 6));
            spacerPanel[i].setBorder(BorderFactory.createMatteBorder(0,0,1,0,Color.DARK_GRAY));
            spacerPanel[i].setOpaque(false);
        }
        JPanel[] dataPanel = new JPanel[lblarr.length];
        for (int i = 0; i < dataPanel.length; i++)
        {
            dataPanel[i] = new JPanel();
            dataPanel[i].setPreferredSize(new Dimension(230, 22));
            dataPanel[i].setOpaque(false);
            lblarr[i].setPreferredSize(new Dimension(30, 18));
            dataPanel[i].add(lblarr[i]);
            txtCarr[i][0].setPreferredSize(new Dimension(90, 18));
            txtCarr[i][1].setPreferredSize(new Dimension(90, 18));
            dataPanel[i].add(txtCarr[i][0]);
            dataPanel[i].add(txtCarr[i][1]);
        }

        JPanel alphaPanel = new JPanel();
        alphaPanel.setOpaque(false);
        alphaPanel.setPreferredSize(new Dimension(145, 24));
        JLabel lblalpha = new JLabel("alpha");
        lblalpha.setPreferredSize(new Dimension(40, 18));
        alphaPanel.add(lblalpha);
        txtalpha.setPreferredSize(new Dimension(90, 18));
        txtalpha.setEditable(false);
        alphaPanel.add(txtalpha);

        JPanel thetaPanel = new JPanel();
        thetaPanel.setOpaque(false);
        thetaPanel.setPreferredSize(new Dimension(145, 24));
        JLabel lbltheta = new JLabel("theta");
        lbltheta.setPreferredSize(new Dimension(40, 18));
        thetaPanel.add(lbltheta);
        txttheta.setPreferredSize(new Dimension(90, 18));
        txttheta.setEditable(false);
        thetaPanel.add(txttheta);

        JPanel aPanel = new JPanel();
        aPanel.setOpaque(false);
        aPanel.setPreferredSize(new Dimension(145, 24));
        JLabel lbla = new JLabel("a");
        lbla.setPreferredSize(new Dimension(40, 18));
        aPanel.add(lbla);
        txta.setPreferredSize(new Dimension(90, 18));
        aPanel.add(txta);

        JPanel bPanel = new JPanel();
        bPanel.setOpaque(false);
        bPanel.setPreferredSize(new Dimension(145, 24));
        JLabel lblb = new JLabel("b");
        lblb.setPreferredSize(new Dimension(40, 18));
        bPanel.add(lblb);
        txtb.setPreferredSize(new Dimension(90, 18));
        bPanel.add(txtb);

        JPanel startPanel = new JPanel();
        startPanel.setOpaque(false);
        startPanel.setPreferredSize(new Dimension(185, 24));
        JLabel lblstart = new JLabel("start");
        lblstart.setPreferredSize(new Dimension(40, 18));
        startPanel.add(lblstart);
        txtstartx.setPreferredSize(new Dimension(60, 18));
        startPanel.add(txtstartx);
        txtstarty.setPreferredSize(new Dimension(60, 18));
        startPanel.add(txtstarty);

        JPanel rangePanel = new JPanel();
        rangePanel.setOpaque(false);
        rangePanel.setPreferredSize(new Dimension(185, 24));
        JLabel lblrange = new JLabel("radius");
        lblrange.setPreferredSize(new Dimension(50, 18));
        rangePanel.add(lblrange);
        txtrange.setPreferredSize(new Dimension(70, 18));
        rangePanel.add(txtrange);

        JPanel printPanel = new JPanel();
        printPanel.setOpaque(false);
        printPanel.setPreferredSize(new Dimension(185, 24));
        invertChk.setOpaque(false);
        printPanel.add(invertChk);
        printChk.setOpaque(false);
        printPanel.add(printChk);

        JPanel overridePanel = new JPanel();        // temporary special panel just for the purpose of overriding
        overridePanel.setOpaque(false);             // beta and PNG Color
        overridePanel.setPreferredSize(new Dimension(185, 24));
        txtbeta_override.setPreferredSize(new Dimension(80, 18));
        overridePanel.add(txtbeta_override);
        txtcolor_override.setPreferredSize(new Dimension(80, 18));
        overridePanel.add(txtcolor_override);

        final JPanel parmsPanel = new JPanel();
        parmsPanel.setBackground(new Color(200, 221, 242));
        parmsPanel.add(alphaPanel);
        parmsPanel.add(thetaPanel);
        parmsPanel.add(aPanel);
        parmsPanel.add(bPanel);
        parmsPanel.add(spacerPanel[0]);
        for (int i = 0; i < dataPanel.length; i++)
            parmsPanel.add(dataPanel[i]);
        parmsPanel.add(spacerPanel[1]);
        parmsPanel.add(btnCalc);
        parmsPanel.add(btnResetHenon);
        parmsPanel.add(startPanel);
        parmsPanel.add(rangePanel);
        parmsPanel.add(btnReset);
        parmsPanel.add(btnClear);
        parmsPanel.add(btnInc);
        parmsPanel.add(printPanel);
        parmsPanel.add(overridePanel);                          // temporary override of beta and Color
        parmsPanel.setMaximumSize(new Dimension(250, 3000));
        parmsPanel.setPreferredSize(new Dimension(250, 3000));
        parmsPanel.setBorder(BorderFactory.createEtchedBorder());

        final JPanel scatterPanel = new JPanel();
        scatterPanel.setOpaque(false);
        scatterPanel.add(lblImage);

        load_curve_fit();
        init();
        x0 = Double.parseDouble(txtstartx.getText());
        y0 = Double.parseDouble(txtstarty.getText());
        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.X_AXIS));
        getContentPane().setBackground(new Color(200, 221, 242));
        getContentPane().add(parmsPanel);
        getContentPane().add(scatterPanel);
        lblImage.addMouseListener(new MouseAdapter()
        {
            @Override public void mousePressed(MouseEvent e)
            {
                //System.out.println(e.getX() + ", " + e.getY() + ", " + e.getButton() + ", " + e.getModifiers());
                x0 = (2.0*e.getX()/image.getWidth() - 1.0)*Double.parseDouble(txtrange.getText());
                y0 = (1.0 - 2.0*e.getY()/image.getHeight())*Double.parseDouble(txtrange.getText());
                //y0 = x0*Math.tan(H_alpha/2);    // override: snap onto plane of symmetry
                //x0 = 0.4641;
                //y0 = 0.775*x0;                                  // fix fix force slope
                System.out.println("Mouse Clicked (" + x0 + ", " + y0 + ")");
                //System.out.println("Mouse Clicked (" + xclick + ", " + yclick + ")");
                if (e.getX() > 0 && e.getX() < image.getWidth() && e.getY() > 0 && e.getY() < image.getHeight())
                {
                    //image.setRGB(e.getX(), e.getY(), Color.BLUE.getRGB());
                    draw_cross (x0, y0, Color.LIGHT_GRAY);
                    lblImage.repaint();
                }
            }
        });

        setVisible(true);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);

        btnCalc.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                btnCalc.setEnabled(false);
                for (int i = 0; i < txtCarr.length; i++)
                    for (int j = 0; j < 2; j++)
                        if (txtCarr[i][j].getText().isEmpty())
                            txtCarr[i][j].setText("0");
                refresh_graph();
                btnCalc.setEnabled(true);
            }
        });

        btnResetHenon.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                // overwrite Henon 'alpha' and 'beta'   TEST CODE ONLY, DO NOT KEEP!!!
                //Henon_a = Double.parseDouble(txtCarr[0][0].getText());
                //Henon_b = Double.parseDouble(txtCarr[0][1].getText());
                //System.out.println("Henon a b =, " + Henon_a + ", " + Henon_b);
            }
        });

        btnReset.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                x0 = Double.parseDouble(txtstartx.getText());
                y0 = Double.parseDouble(txtstarty.getText());
                //y0 = x0*Math.tan(H_alpha/2);    // override: snap onto plane of symmetry
                //y0 = x0*0.23/0.674;                 // aim for sixfold island resonance at cos = .40
                //txtstarty.setText(String.format("%.6f", y0));       // truncate display ONLY
                //System.out.println("btnReset ," + Delay_a + ", " + x0 + ", " + y0);
                //x0 = xclick;
                //y0 = yclick;
                Henon_b = Double.parseDouble(txtbeta_override.getText());                // temporary override
                System.out.println("over-ride, " + x0 + ", " + y0 + ", " + Henon_b + ", " + txtcolor_override.getText());
            }
        });

        btnClear.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                DC.clearRect(0, 0, image.getWidth(), image.getHeight());
                lblImage.repaint();
            }
        });
        //transform_circle();
        //System.out.println("H_alpha/2 = " + H_alpha/2*180/Math.PI);
        btnInc.addActionListener(new AbstractAction()
        {
            public void actionPerformed(ActionEvent event)
            {
                //Delay_a = 0.0001*Math.round(10000*Delay_a - 1);
                //Delay_a = 0.00001*Math.round(100000*Delay_a - 1);
                //System.out.println("Delay_a = " + Delay_a);
                Henon_b = 0.00001*Math.round(100000*Henon_b + 1);
                System.out.println("Henon_b = " + Henon_b);
            }
        });
        txtrange.addKeyListener(new KeyAdapter()
        {
            @Override public void keyPressed(KeyEvent e)
            {
                if (e.getKeyCode() == KeyEvent.VK_ESCAPE)
                {
                    String fname = "Henon_scatter_" + Henon_a + "_" + Henon_b + "_" + Henon_R + "_" + Henon_S;
                    Main.save_PNG(image, fname);
                }
            }
        });
    }

    protected static void init()
    {
        System.out.println("\ncall init from Chua_Simul_3");
        Chua_iTstr  = hdr.split(",")[0].trim();
        Chua_alpha  = Double.parseDouble(hdr.split(",")[1]);
        Chua_beta   = Double.parseDouble(hdr.split(",")[2]);
        Chua_gamma  = Double.parseDouble(hdr.split(",")[3]);
        Chua_a      = Double.parseDouble(hdr.split(",")[4]);
        Chua_c      = Double.parseDouble(hdr.split(",")[5]);
        Chua_Period = Double.parseDouble(hdr.split(",")[6]);
        Chua_delt   = Double.parseDouble(hdr.split(",")[7]);
        //for (int i = 0; i < Cx.length; i++)
        //    Cx[i] = Double.parseDouble(hdr.split(",")[8 + i]);
        //for (int i = 0; i < Cy.length; i++)
        //    Cy[i] = Double.parseDouble(hdr.split(",")[18 + i]);
        for (int i = 0; i < C_g_2_R.length; i++)
            for (int j = 0; j < C_g_2_R.length; j++)
            {
                C_g_2_R[i][j] /= 2.0;
                C_g_2_I[i][j] /= 2.0;
            }
        for (int i = 0; i < C_g_3_R.length; i++)
            for (int j = 0; j < C_g_3_R.length; j++)
            {
                C_g_3_R[i][j] /= 4.0;
                C_g_3_I[i][j] /= 4.0;
            }

        // define Taylor coeff 'g' using Kuznetsov notation and a 2D point for complex numbers

        Point2D.Double mu10 = new Point2D.Double(Double.parseDouble(hdr.split(",")[9]), -Double.parseDouble(hdr.split(",")[10])); // first-order response
        //System.out.println("linear a b, " + mu10);
        Point2D.Double mu01 = Main.conjugate(mu10);                  // conjugate of mu
        Point2D.Double mu20 = Main.multiply(mu10, mu10);
        Point2D.Double mu11 = Main.multiply(mu10, mu01);
        Point2D.Double mu02 = Main.multiply(mu01, mu01);
        Point2D.Double mu30 = Main.multiply(mu20, mu10);
        Point2D.Double mu21 = Main.multiply(mu20, mu01);
        Point2D.Double mu12 = Main.multiply(mu11, mu01);
        Point2D.Double mu03 = Main.multiply(mu02, mu01);
        Point2D.Double[] g_2 = new Point2D.Double[3];           // quadratic response
        Point2D.Double[] g_3 = new Point2D.Double[4];           // cubic response
        //System.out.println("a b, " + mu10);
        for (int i = 0; i < g_2.length; i++)
        {
            g_2[i] = new Point2D.Double(0.0, 0.0);
            for (int j = 0; j < g_2.length; j++)
            {
                g_2[i].x += C_g_2_R[i][j]*Carr[j][0] - C_g_2_I[i][j]*Carr[j][1];
                g_2[i].y += C_g_2_R[i][j]*Carr[j][1] + C_g_2_I[i][j]*Carr[j][0];
            }
        }
        for (int i = 0; i < g_3.length; i++)
        {
            g_3[i] = new Point2D.Double(0.0, 0.0);
            for (int j = 0; j < g_3.length; j++)
            {
                g_3[i].x += C_g_3_R[i][j]*Carr[j + 3][0] - C_g_3_I[i][j]*Carr[j + 3][1];
                g_3[i].y += C_g_3_R[i][j]*Carr[j + 3][1] + C_g_3_I[i][j]*Carr[j + 3][0];
            }
        }
        //System.out.println("convert_Cxy_to_gij:" + hdr);
        System.out.println("convert_Cxy_to_gij: " + Chua_iTstr + ", " + Chua_alpha + ", " + Chua_beta + ", " + Chua_gamma + ", " + Chua_a + ", " + Chua_c + ", " + Chua_Period + ", " + Chua_delt);
        //System.out.print("Cx ");
        //for (int i = 0; i < Cx.length; i++)
        //    System.out.print(", " + Cx[i]);
        //System.out.println();
        //System.out.print("Cy ");
        //for (int i = 0; i < Cy.length; i++)
        //    System.out.print(", " + Cy[i]);
        //System.out.println();
        //double xinv = Carr[2][0]/(Carr[1][0]*Carr[1][0] - 4*Carr[0][0]*Carr[2][0]);  // extremum x of the inverse transform
        //double yinv = Carr[0][1]/(Carr[1][1]*Carr[1][1] - 4*Carr[0][1]*Carr[2][1]);  // extremum x of the inverse transform
        //System.out.println("extremum of the inverse , " + xinv + ", " + yinv + ", " + (Carr[1][0]*Carr[1][0] - 4*Carr[0][0]*Carr[2][0]) + ", " + (Carr[1][1]*Carr[1][1] - 4*Carr[0][1]*Carr[2][1]));

        System.out.print("g_2");
        for (int i = 0; i < g_2.length; i++)
            System.out.print(", " + g_2[i]);
        System.out.println();
        System.out.print("g_3");
        for (int i = 0; i < g_3.length; i++)
            System.out.print(", " + g_3[i]);
        System.out.println();

        // coeff 'h' of nonlinear transform

        Point2D.Double htemp = new Point2D.Double();            // (Kuznetsov p. 151)
        htemp = Main.subtract(mu20, mu10);
        h_2[0] = Main.divide(g_2[0], htemp);
        htemp = Main.subtract(mu11, mu10);
        h_2[1] = Main.divide(g_2[1], htemp);
        htemp = Main.subtract(mu02, mu10);
        h_2[2] = Main.divide(g_2[2], htemp);
        System.out.print("h_2");
        for (int i = 0; i < h_2.length; i++)
            System.out.print(", " + h_2[i]);
        System.out.println();

        Point2D.Double[] h_3 = new Point2D.Double[4];           // cubic response
        htemp = Main.subtract(mu30, mu10);                      // (Kuznetsov p. 152)
        h_3[0] = Main.divide(g_3[0], htemp);
        htemp = Main.subtract(mu21, mu10);
        h_3[1] = Main.divide(g_3[1], htemp);
        htemp = Main.subtract(mu12, mu10);
        h_3[2] = Main.divide(g_3[2], htemp);
        htemp = Main.subtract(mu03, mu10);
        h_3[3] = Main.divide(g_3[3], htemp);
        System.out.print("h_3");
        for (int i = 0; i < h_3.length; i++)
            System.out.print(", " + h_3[i]);
        System.out.println();

        Point2D.Double htemp2 = new Point2D.Double();           // terms A, B, C, D
        Point2D.Double htemp3 = new Point2D.Double();           // N_S bifurc, Dec 21, 2021
        Point2D.Double htemp4 = new Point2D.Double();
        Point2D.Double htemp5 = new Point2D.Double();
        htemp  = Main.multiply(0.5, Main.multiply(h_2[0], h_2[0]));
        htemp2 = Main.multiply(0.5, Main.multiply(h_2[1], Main.conjugate(h_2[2])));
        htemp3 = Main.multiply(-1.0/6.0, h_3[0]);                   // add cubic term h30
        //htemp3 = new Point2D.Double(0, 0);                        // DISABLE cubic TEMPORARY
        inv_3[0] = Main.add(htemp, Main.add(htemp2, htemp3));             // term A

        htemp  = Main.multiply(1.5, Main.multiply(h_2[0], h_2[1]));
        htemp2 = Main.multiply(0.5, Main.multiply(h_2[2], Main.conjugate(h_2[2])));
        htemp3 = Main.multiply(h_2[1], Main.conjugate(h_2[1]));     // do NOT add cubic term h21
        inv_3[1] = Main.add(htemp, Main.add(htemp2, htemp3));             // term B

        htemp  = Main.multiply(0.5, Main.multiply(h_2[0], h_2[2]));
        htemp2 = Main.multiply(h_2[1], h_2[1]);
        htemp3 = Main.multiply(0.5, Main.multiply(h_2[1], Main.conjugate(h_2[0])));
        htemp4 = Main.multiply(h_2[2], Main.conjugate(h_2[1]));
        htemp5 = Main.multiply(-1.0/2.0, h_3[2]);                   // add cubic term h12
        //htemp5 = new Point2D.Double(0, 0);                        // DISABLE cubic TEMPORARY
        inv_3[2] = Main.add(htemp, Main.add(htemp2, Main.add(htemp3, Main.add(htemp4, htemp5))));    // term C

        htemp  = Main.multiply(0.5, Main.multiply(h_2[1], h_2[2]));
        htemp2 = Main.multiply(0.5, Main.multiply(h_2[2], Main.conjugate(h_2[0])));
        htemp3 = Main.multiply(-1.0/6.0, h_3[3]);                   // add cubic term h03
        //htemp3 = new Point2D.Double(0, 0);                        // DISABLE cubic TEMPORARY
        inv_3[3] = Main.add(htemp, Main.add(htemp2, htemp3));             // term D

        System.out.print("inv_3");
        for (int i = 0; i < inv_3.length; i++)
            System.out.print(", " + inv_3[i]);
        System.out.println();
        //System.out.println("test transform : " + transform(new Point2D.Double(-0.2374, 1.825)));
    }

//    protected static Point2D.Double simulate_cubic(double x, double y)
//    {
//        // use Cx, Cy from hdr to simulate a cubic map
//        //return new Point2D.Double(Cx[1]*x + Cx[2]*y
//        //                        , Cy[1]*x + Cy[2]*y);
//        return new Point2D.Double(Cx[1]*x + Cx[2]*y + Cx[3]*x*x + Cx[4]*x*y + Cx[5]*y*y + Cx[6]*x*x*x + Cx[7]*x*x*y + Cx[8]*x*y*y + Cx[9]*y*y*y
//                                , Cy[1]*x + Cy[2]*y + Cy[3]*x*x + Cy[4]*x*y + Cy[5]*y*y + Cy[6]*x*x*x + Cy[7]*x*x*y + Cy[8]*x*y*y + Cy[9]*y*y*y);
//    }

    private void transform_circle()
    {
        double x, y;
        Point2D.Double trans;
        for (int i = 0; i < 60; i++)
        {
            x = Math.cos(6*i*Math.PI/180);
            y = Math.sin(6*i*Math.PI/180);
            trans = transform(x, y);
            System.out.println(6*i + ", " + x + ", " + y + ", " + trans.x + ", " + trans.y);
        }
    }

    private static void draw_cross(double xc, double yc, Color clr)
    {
        DC.setColor(clr);
        int setX = (int) ((1 + xc/Double.parseDouble(txtrange.getText()))*image.getWidth()/2);
        int setY = (int) ((1 - yc/Double.parseDouble(txtrange.getText()))*image.getWidth()/2);
        DC.drawLine(setX - 2, setY, setX + 2, setY);
        DC.drawLine(setX, setY - 2, setX, setY + 2);
    }

    protected static Point2D.Double transform(double x, double y)
    {
        // transform z to w
        Point2D.Double z10 = new Point2D.Double(x, y);
        Point2D.Double z01 = Main.conjugate(z10);
        Point2D.Double z20 = Main.multiply(z10, z10);
        Point2D.Double z11 = Main.multiply(z10, z01);
        Point2D.Double z02 = Main.multiply(z01, z01);
        Point2D.Double z30 = Main.multiply(z20, z10);
        Point2D.Double z21 = Main.multiply(z20, z01);
        Point2D.Double z12 = Main.multiply(z11, z01);
        Point2D.Double z03 = Main.multiply(z02, z01);

        z10 = Main.subtract(z10, Main.multiply(0.5, h_2[0], z20));
        z10 = Main.subtract(z10, Main.multiply(h_2[1], z11));
        z10 = Main.subtract(z10, Main.multiply(0.5, h_2[2], z02));
        z10 = Main.add(z10, Main.multiply(inv_3[0], z30));
        z10 = Main.add(z10, Main.multiply(inv_3[1], z21));
        z10 = Main.add(z10, Main.multiply(inv_3[2], z12));
        z10 = Main.add(z10, Main.multiply(inv_3[3], z03));
        return z10;
    }

//    private void plot_Logistic_boundary()
//    {
//        // see: Pounder & Rogers, region E0, p. 560
//        for (int i = 0; i < (1.0 - 1.0/Delay_a)*image.getWidth()/2; i++)
//            //image.setRGB(image.getWidth()/2  + i, image.getWidth()/2  - (int) (1.5), Color.LIGHT_GRAY.getRGB());
//            image.setRGB(image.getWidth()/2  + i, image.getWidth()/2  - (int) (image.getWidth()/2.0/Delay_a/(1.0 - 2.0*i/image.getWidth())), Color.LIGHT_GRAY.getRGB());
//    }

    protected static Point2D.Double rotate(Point2D.Double pt_in, double rot)
    {
        // rotate point by angle : Henon_1969 Eq. 13
        return new Point2D.Double(pt_in.x*Math.cos(rot) - pt_in.y*Math.sin(rot),
                                  pt_in.x*Math.sin(rot) + pt_in.y*Math.cos(rot));
    }

    protected static Point2D.Double skew(Point2D.Double pt_in, double skew)
    {
        // skew point: y = y - skew*x*x : Henon_1969 Eq. 12
        return new Point2D.Double(pt_in.x, pt_in.y - skew*pt_in.x*pt_in.x);
    }

    protected static Point2D.Double reflect(Point2D.Double pt_in, double alpha)
    {
        // reflect point about a line at angle alpha (Henon_1969)
        return new Point2D.Double(pt_in.x*Math.cos(2*alpha) + pt_in.y*Math.sin(2*alpha), pt_in.x*Math.sin(2*alpha) - pt_in.y*Math.cos(2*alpha));
    }
/*
    protected static Point2D.Double map_transform()
    {
        // implement Henon map : Henon_1969 Eq. 4-8
        // scale C and E according to D in Eq.7
        //double beta = -H_alpha/2 + Math.PI/2;       // even in x
        double beta = -H_alpha + Math.PI;               // original x^2 model

        double D = -2*Math.sin(H_alpha + beta)*Math.cos(H_alpha + beta);
        double C = Math.cos(H_alpha + beta)*Math.cos(H_alpha + beta);
        double E = Math.sin(H_alpha + beta)*Math.sin(H_alpha + beta);
        double R = C*x0*x0 + D*x0*y0 + E*y0*y0;
        //System.out.println("C D E R, " + C + ", " + D + ", " + E + ", " + R);
        return new Point2D.Double(x0*Math.cos(H_alpha) - y0*Math.sin(H_alpha) + R*Math.sin(beta),
                                  x0*Math.sin(H_alpha) + y0*Math.cos(H_alpha) + R*Math.cos(beta));
    }
*/
    private void refresh_graph()
    {
        PrintWriter fout = null;
        double a = Double.parseDouble(txta.getText());
        double b = Double.parseDouble(txtb.getText());
        double range = Double.parseDouble(txtrange.getText());
        Point2D.Double ztrans;
        //Point2D.Double pt;
        double xplt, yplt;              // transformed coord
        //double H_del = 1*(Math.PI/2 - H_alpha/2);         // rotate Henon by angle del
        double x, y;
        //double x2 = 0.6575379230230298;                              // initialize x at t-2
        int N = 30000;
        //int N = 29982;                  //

        //x0 = 0.46040285344901516;
        //x0 = -0.333; //0.5076983932892956;                               // FIX FIX temporary override
        //y0 = -0.275; //-0.43;
        //pt = new Point2D.Double(-0.333, -0.275);
        //pt = new Point2D.Double(-0.2760382492972518, -0.3180380523020866);
        //x0 = pt.x;
        //y0 = pt.y;
        //x0 = 0.15; // 0.28;                  // FIX FIX temporary override
        //y0 = 0.11625; // 0.46;
        //H_alpha = Math.acos(H_cos);

        System.out.println("refresh_graph, " + iT + ", " + invertChk.isSelected() + ", " + x0 + ", " + y0);
        //System.out.println("refresh_graph, " + iT + ", " + Delay_a + ", " + x0 + ", " + y0);
        for (int i = 0; i < txtCarr.length; i++)
        {
            Carr[i][0] = Double.parseDouble(txtCarr[i][0].getText());   // Cx
            Carr[i][1] = Double.parseDouble(txtCarr[i][1].getText());   // Cy
        }
        if (printChk.isSelected())
            try
            {
                //String fname = "Chua_Simul_scatter_" + hdr.split(",")[0].trim() + "_" + hdr.split(",")[1].trim() + "_" + (hdr.split(",").length/2 - 7);
                String fname = "Henon_scatter_" + Henon_a + "_" + Henon_b + "_" + Henon_R + "_" + Henon_S;
                //String fname = "Henon_1969_scatter_" + H_cos + "_" + String.format("%.6f", x0) + "_" + String.format("%.6f", y0);
                //String fname = "Henon_theory_no4_g50_g32_g14_" + hdr.split(",")[0].trim() + "_" + hdr.split(",")[1].trim() + "_" + hdr.split(",")[2].trim() + "_" + hdr.split(",")[3].trim() + "_" + (hdr.split(",").length/2 - 7);
                //String fname = "Henon_quadratic_" + hdr.split(",")[0].trim() + "_" + hdr.split(",")[1].trim() + "_" + hdr.split(",")[2].trim() + "_" + hdr.split(",")[3].trim() + "_" + (hdr.split(",").length/2 - 7);
                //String fname = "DL_Linearized_" + hdr.split(",")[1].trim();
                //String fname = "DL_final_g21_no_quartic_no_quintic_" + hdr.split(",")[1].trim();
                //String fname = "DL_raw_" + Delay_a;
                boolean fexists = new File("C:\\Windows\\Temp\\" + fname + ".csv").exists();
                System.out.println("fexists = " + fexists + ", " + iT);
                FileWriter fw = new FileWriter("C:\\Windows\\Temp\\" + fname + ".csv", true);
                fout = new PrintWriter(fw);
                if (!fexists)
                {
                    fout.println("Cxi, " + Carr[0][0] + ", " + Carr[1][0] + ", " + Carr[2][0] + ", " + Carr[3][0] + ", " + Carr[4][0] + ", " + Carr[5][0] + ", " + Carr[6][0]
                                  + ", " + Carr[0][1] + ", " + Carr[1][1] + ", " + Carr[2][1] + ", " + Carr[3][1] + ", " + Carr[4][1] + ", " + Carr[5][1] + ", " + Carr[6][1]);
                    //fout.println("simul3, " + Chua_alpha + ", " + Chua_beta + ", " + Chua_gamma + ", " + Chua_a + ", " + Chua_c + ", " + Chua_delt);
                    //fout.println("raw   , " + Delay_a + ", NaN, NaN, NaN, NaN, NaN");
                    fout.println("Henon, " + Henon_a + ", " + Henon_b + ", " + Henon_R + ", " + Henon_S + ", " + "NaN" + ", " + "NaN");
                    //fout.println("Henon 1969, " + H_cos + ", " + String.format("%.8f", x0) + ", " + String.format("%.8f", y0) + ", " + String.format("%.8f", 2*Math.tan(H_alpha/2)) + ", " + String.format("%.8f", 2*Math.tan(H_alpha/2)*Math.tan(H_alpha/2)) + ", " + "NaN");
                    //fout.println("Henon 1969, " + H_cos + ", " + String.format("%.8f", x0) + ", " + String.format("%.8f", y0) + ", " + "NaN" + ", " + "NaN" + ", " + "NaN");
                    fout.println("scatter hdr, " + hdr.split(",")[0] + ", " + Math.sqrt(a*a + b*b) + ", " + txttheta.getText());
                    fout.println("init x0 y0, 0, 0, NaN, NaN");         // normal simulation
                    //fout.println("init x0 y0, 0.5, 0.5, NaN, NaN");   // Delayed Logistic
                    fout.println("iter       , x', y'");
                    fout.println(iT + ", " + x0 + ", " + y0);
                }
            }
            catch (java.io.IOException e)
                {System.out.println("Chua_Output.csv save error = " + e);}
        //System.out.println("Carr len = " + Carr.length);
        //double accumx = 0;              // accum x to check for cyclic behaviour
        for (int i = 0; i < N; i++)
        {
            x = a*x0 - b*y0
              + Carr[0][0]*x0*x0 + Carr[1][0]*x0*y0 + Carr[2][0]*y0*y0
              + Carr[3][0]*x0*x0*x0 + Carr[4][0]*x0*x0*y0 + Carr[5][0]*x0*y0*y0 + Carr[6][0]*y0*y0*y0;
            if (Carr.length > 7)        // 7 = cubic
                x = x + Carr[7][0]*x0*x0*x0*x0 + Carr[8][0]*x0*x0*x0*y0 + Carr[9][0]*x0*x0*y0*y0 + Carr[10][0]*x0*y0*y0*y0 + Carr[11][0]*y0*y0*y0*y0;
            if (Carr.length > 12)       // 12 = quartic
                x = x + Carr[12][0]*x0*x0*x0*x0*x0 + Carr[13][0]*x0*x0*x0*x0*y0 + Carr[14][0]*x0*x0*x0*y0*y0 + Carr[15][0]*x0*x0*y0*y0*y0 + Carr[16][0]*x0*y0*y0*y0*y0 + Carr[17][0]*y0*y0*y0*y0*y0;
            if (Carr.length > 18)       // 18 = quintic
                x = x + Carr[18][0]*x0*x0*x0*x0*x0*x0 + Carr[19][0]*x0*x0*x0*x0*x0*y0 + Carr[20][0]*x0*x0*x0*x0*y0*y0 + Carr[21][0]*x0*x0*x0*y0*y0*y0 + Carr[22][0]*x0*x0*y0*y0*y0*y0 + Carr[23][0]*x0*y0*y0*y0*y0*y0 + Carr[24][0]*y0*y0*y0*y0*y0*y0;
            if (Carr.length > 25)
                x = x + Carr[25][0]*x0*x0*x0*x0*x0*x0*x0 + Carr[26][0]*x0*x0*x0*x0*x0*x0*y0 + Carr[27][0]*x0*x0*x0*x0*x0*y0*y0 + Carr[28][0]*x0*x0*x0*x0*y0*y0*y0 + Carr[29][0]*x0*x0*x0*y0*y0*y0*y0 + Carr[30][0]*x0*x0*y0*y0*y0*y0*y0 + Carr[31][0]*x0*y0*y0*y0*y0*y0*y0 + Carr[32][0]*y0*y0*y0*y0*y0*y0*y0;
            y = b*x0 + a*y0
              + Carr[0][1]*x0*x0 + Carr[1][1]*x0*y0 + Carr[2][1]*y0*y0
              + Carr[3][1]*x0*x0*x0 + Carr[4][1]*x0*x0*y0 + Carr[5][1]*x0*y0*y0 + Carr[6][1]*y0*y0*y0;
            if (Carr.length > 7)        // 7 = cubic
                y = y + Carr[7][1]*x0*x0*x0*x0 + Carr[8][1]*x0*x0*x0*y0 + Carr[9][1]*x0*x0*y0*y0 + Carr[10][1]*x0*y0*y0*y0 + Carr[11][1]*y0*y0*y0*y0;
            if (Carr.length > 12)
                y = y + Carr[12][1]*x0*x0*x0*x0*x0 + Carr[13][1]*x0*x0*x0*x0*y0 + Carr[14][1]*x0*x0*x0*y0*y0 + Carr[15][1]*x0*x0*y0*y0*y0 + Carr[16][1]*x0*y0*y0*y0*y0 + Carr[17][1]*y0*y0*y0*y0*y0;
            if (Carr.length > 18)
                y = y + Carr[18][1]*x0*x0*x0*x0*x0*x0 + Carr[19][1]*x0*x0*x0*x0*x0*y0 + Carr[20][1]*x0*x0*x0*x0*y0*y0 + Carr[21][1]*x0*x0*x0*y0*y0*y0 + Carr[22][1]*x0*x0*y0*y0*y0*y0 + Carr[23][1]*x0*y0*y0*y0*y0*y0 + Carr[24][1]*y0*y0*y0*y0*y0*y0;
            if (Carr.length > 25)
                y = y + Carr[25][1]*x0*x0*x0*x0*x0*x0*x0 + Carr[26][1]*x0*x0*x0*x0*x0*x0*y0 + Carr[27][1]*x0*x0*x0*x0*x0*y0*y0 + Carr[28][1]*x0*x0*x0*x0*y0*y0*y0 + Carr[29][1]*x0*x0*x0*y0*y0*y0*y0 + Carr[30][1]*x0*x0*y0*y0*y0*y0*y0 + Carr[31][1]*x0*y0*y0*y0*y0*y0*y0 + Carr[32][1]*y0*y0*y0*y0*y0*y0*y0;

            // TEMPORARY over-ride with a generalized Henon map (Gonchenko)
            if (i < 5)
                System.out.println(i + ", " + Henon_a + ", " + Henon_b + ", " + x0 + ", " + y0);
            x = y0;
            y = Henon_a - Henon_b*x0 - y0*y0 + Henon_R*x0*y0 + Henon_S*y0*y0*y0;
            // END of over-ride

            // TEMPORARY over-ride with inverse of generalized Henon map (Gonchenko)
            //y = x0;
            //x = (Henon_a - x0*x0 - y0)/(Henon_b - Henon_R*x0);
            // END of over-ride

            // TEMPORARY over-ride with a Henon map (Arrowsmith p. 54) : Henon Eq.11
            //x = x0*Math.cos(H_alpha) - y0*Math.sin(H_alpha) + x0*x0*Math.sin(H_alpha);    // org code x'*x'
            //y = x0*Math.sin(H_alpha) + y0*Math.cos(H_alpha) - x0*x0*Math.cos(H_alpha);
            //x = x0*Math.cos(H_alpha) - y0*Math.sin(H_alpha) + y0*y0*Math.cos(H_alpha);      // org code y'*y'
            //y = x0*Math.sin(H_alpha) + y0*Math.cos(H_alpha) + y0*y0*Math.sin(H_alpha);

            // TEMPORARY over-ride with a Henon map (Arrowsmith p. 54) : Henon Eq.11
            // general map as per p.17 of insert 'Henon_1969_map' Feb 20, 2026
            //double H_beta = -H_alpha - Math.PI/4;           // beta is an adjustable rotation
            //double H_beta = -H_alpha - Math.PI/2;
            //double R = x0*Math.cos(H_alpha + H_beta) - y0*Math.sin(H_alpha + H_beta);
            //x = x0*Math.cos(H_alpha) - y0*Math.sin(H_alpha) - R*R*Math.sin(H_beta);     // general rotation angle
            //y = x0*Math.sin(H_alpha) + y0*Math.cos(H_alpha) - R*R*Math.cos(H_beta);

            //x = x0*Math.cos(H_alpha) - y0*Math.sin(H_alpha) + y0*y0*Math.cos(H_alpha);  // test code y'*y'
            //y = x0*Math.sin(H_alpha) + y0*Math.cos(H_alpha) + y0*y0*Math.sin(H_alpha);
            //double H_beta = -0.99*H_alpha;
            //double R = (1 - Math.cos(4*(H_alpha + H_beta)))*x0*x0/2 + y0*y0;
            //x = x0*Math.cos(H_alpha) - y0*Math.sin(H_alpha) + Math.cos(H_alpha)*R; // new code y'*y'
            //y = x0*Math.sin(H_alpha) + y0*Math.cos(H_alpha) + Math.sin(H_alpha)*R;
            // END of over-ride

            //pt = map_transform();
            //x = pt.x;
            //y = pt.y;

            // over-ride with a Henon map (Arrowsmith) (rotated CCW by (90 - alpha/2)
            //double x_trans = x0*Math.cos(H_del) + y0*Math.sin(H_del);       // back-transform x'
            //x = x0*Math.cos(H_alpha) - y0*Math.sin(H_alpha) + x_trans*x_trans*Math.sin(H_alpha + H_del);
            //y = x0*Math.sin(H_alpha) + y0*Math.cos(H_alpha) - x_trans*x_trans*Math.cos(H_alpha + H_del);
            // END of over-ride

            // TEMPORARY over-ride with the inverse of a Henon map (Arrowsmith p. 54) : Henon Eq.14
            //double x_trans =  x0*Math.cos(H_alpha - H_del) + y0*Math.sin(H_alpha - H_del);
            //x =  x0*Math.cos(H_alpha) + y0*Math.sin(H_alpha) + Math.sin(H_del)*x_trans*x_trans;
            //y = -x0*Math.sin(H_alpha) + y0*Math.cos(H_alpha) + Math.cos(H_del)*x_trans*x_trans;
            // END of over-ride

            // over-ride with a Henon map (Arrowsmith) (rotated CCW by (90 - alpha/2) (DO NOT USE - WRONG)
            // expressed as a second-order d.e. in x, in the time sequence : (x2, x0, x)
            // see: "Fourier Analysis of Henon 'sin_cos' map, p.13"
            //x = 2*x0*Math.cos(H_alpha) - x2 + Math.cos(H_alpha/2)*(x0*x0 - x2*x2);
            //y = (-x*Math.cos(H_alpha) + x0 + x2*x2*Math.sin(Math.PI/2 - H_alpha/2))/Math.sin(H_alpha);

            // second-order d.e. in x, as above, but use formulas from p.16 (decouple Henon_1969)
            //x = 2*x0*Math.cos(H_alpha) - x2 + Math.sin(H_alpha)*x0*x0;
            //y = (-x*Math.cos(H_alpha) + x0)/Math.sin(H_alpha);
            // this algorithm does NOT use y0

            //pt = reflect(pt, H_alpha/2);
            //pt = rotate(skew(reflect(pt, H_alpha/2), 1), H_alpha);       // reflect_skew_rotate Henon map
            //pt = reflect(rotate(skew(reflect(pt, H_alpha/2), 1), H_alpha), H_alpha/2);       // reflect_skew_rotate Henon map
            //pt = rotate(skew(pt, 1), H_alpha);       // implement normal Henon map
            //pt = skew(rotate(pt, -H_alpha), -1);       // implement inverse Henon map
            //pt = rotate(skew(pt, 0.5), H_alpha/2);       // implement half Henon map
            //x = pt.x;
            //y = pt.y;

            // TEMPORARY over-ride with a Delayed Logistic Map (Aronson)
            //x = y0;
            //y = Delay_a*y0*(1 - x0);
            //accumx += x;
            // END of over-ride

            // TEMPORARY over-ride with inverse of a Delayed Logistic Map (Aronson)
            //if (i < 20)
            //    System.out.println(i + ", " + Delay_a + ", " + x0 + ", " + y0);
            //y = x0;
            //x = 1 - y0/Delay_a/x0;
            // END of over-ride

            if (invertChk.isSelected())
            {
                ztrans = transform(x, y);
                xplt = ztrans.x;
                yplt = ztrans.y;
            }
            else
            {
                xplt = x;
                yplt = y;
            }
            if (image.getWidth()/2 + xplt/range*image.getWidth()/2 < 0
            ||  image.getWidth()/2 + xplt/range*image.getWidth()/2 > image.getWidth()
            ||  image.getHeight()/2 - yplt/range*image.getHeight()/2 < 0
            ||  image.getHeight()/2 - yplt/range*image.getHeight()/2 > image.getHeight())
                System.out.println("BAD: " + xplt + ", " + yplt);
            else
                //image.setRGB(image.getWidth()/2 + (int) (xplt/range*image.getWidth()/2), image.getHeight()/2  - (int) (yplt/range*image.getHeight()/2), Color.BLUE.getRGB());
                image.setRGB(image.getWidth()/2 + (int) (xplt/range*image.getWidth()/2), image.getHeight()/2  - (int) (yplt/range*image.getHeight()/2), Integer.parseInt(txtcolor_override.getText()));
            //x2 = x0;                        // use ONLY for the second-order d.e. initilization
            x0 = x;
            y0 = y;
            iT++;
            if (fout != null)
                fout.println(iT + ", " + xplt + ", " + yplt);           // normal output
            //                    + ", " + ( xplt*Math.cos(H_alpha/2 + Math.PI/2) + yplt*Math.sin(H_alpha/2 + Math.PI/2))
            //                    + ", " + (-xplt*Math.sin(H_alpha/2 + Math.PI/2) + yplt*Math.cos(H_alpha/2 + Math.PI/2)));                                 // rotated output
        }
        //System.out.println("accumx ," + Delay_a + ", " + x0 + ", " + accumx);
        //int r = (int) (Math.sqrt(x0*x0 + y0*y0)/range*image.getWidth()/2);
        //DC.setColor(Color.blue);
        //double scale = 2;
        DC.setColor(new Color(192, 128, 96));
        DC.drawLine(0, image.getHeight()/2, image.getWidth(), image.getHeight()/2);
        DC.drawLine(image.getWidth()/2, 0, image.getWidth()/2, image.getHeight());

        DC.setFont(new Font( "SansSerif", Font.BOLD, 12 ));
        DC.setColor(Color.BLACK);                   // = -16777216
        DC.drawString("1.0205", 240 -  75, 240 +  75);
        DC.drawString("1.023", 240 - 130, 240 + 130);
        DC.drawString("1.0249777", 240 - 140, 240 + 175);
        DC.drawString("non-resonant", 360, 20);
        DC.setColor(new Color(255, 128, 0));        // = -32768
        DC.drawString("five sinks", 360, 40);

        //DC.drawOval(image.getHeight()/2 - r, image.getWidth()/2 - r, 2*r, 2*r);
        //DC.drawLine((int) (image.getWidth()/2*(1 - scale*Math.cos(H_alpha/2))), (int) (image.getHeight()/2*(1 + scale*Math.sin(H_alpha/2))),
        //            (int) (image.getWidth()/2*(1 + scale*Math.cos(H_alpha/2))), (int) (image.getHeight()/2*(1 - scale*Math.sin(H_alpha/2))));
        //draw_cross (2*Math.tan(H_alpha/2), 2*Math.tan(H_alpha/2)*Math.tan(H_alpha/2), DC.getColor());
        //image.setRGB(image.getWidth()/2   + (int) (2*Math.tan(H_alpha/2)/range*image.getWidth()/2),
        //             image.getHeight()/2  - (int) (2*Math.tan(H_alpha/2)*Math.tan(H_alpha/2)/range*image.getHeight()/2), Color.MAGENTA.getRGB());
        //plot_Logistic_boundary();
        lblImage.repaint();
        if (fout != null)
            fout.close();
    }

    private static void load_curve_fit()
    {
        //System.out.println("hdr " + hdr.split(",").length);
        double a = Double.parseDouble(hdr.split(",")[9]);
        txta.setText(String.format("%.8f", a));
        double b = -Double.parseDouble(hdr.split(",")[10]);
        txtb.setText(String.format("%.8f", b));
        txtalpha.setText(String.format("%.8f", Math.sqrt(a*a + b*b) - 1));      // first-order response
        txttheta.setText(String.format("%.8f", Math.atan2(b, a)*180/Math.PI));  // (redundant)
        System.out.println("a b, " + a + ", " + b);
        System.out.println();
        Carr = new double[hdr.split(",").length/2 - 7][2];  // use only quadratic and higher coeff
        for (int i = 0; i < Carr.length; i++)
        {
            Carr[i][0] = Double.parseDouble(hdr.split(",")[i + 11]);
            Carr[i][1] = Double.parseDouble(hdr.split(",")[i + 11 + Carr.length + 3]);
            if (i > 24)                     // force a quintic fit (test only !!!)
            {
                Carr[i][0] = 0;
                Carr[i][1] = 0;
            }
            System.out.println("Carr " + i + ", " + Carr[i][0] + ", " + Carr[i][1]);
        }
        for (int i = 0; i < txtCarr.length; i++)
        {
            txtCarr[i][0].setText(String.format("%.6f", Carr[i][0]));
            txtCarr[i][1].setText(String.format("%.6f", Carr[i][1]));
        }
    }

    public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                //Chua_Simul_3.init();
                Chua_Simul_3 dlg = new Chua_Simul_3();
                //scatterPanel.addMouseListener(this);
            }
        });
    }
}
