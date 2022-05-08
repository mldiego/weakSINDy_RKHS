function dx = monoD(t,x) 
dx(1,1) = 0.0006749859334188590409553398785647*x(1)^3*x(2)^2 - 0.00010744270815508938365923086166731*x(1)^2*x(2)^3 - 0.00084138879744855010756054980447516*x(1)^2*x(3)^2 - 0.00012126610372184543074070006696275*x(1)^2*x(2)^2 - 0.00018998784610704810837944478407735*x(1)^3*x(3)^2 + 0.00071817589759093358736663503805175*x(2)^2*x(3)^2 + 0.00021789173663400163150072330608964*x(2)^2*x(3)^3 - 0.0001308905135983806644617288839072*x(1)*x(2)^3 - 0.00025008626734290828608209267258644*x(1)*x(3)^2 - 0.0028846005195801716070036491146311*x(1)*x(3)^3 - 0.00020157475507179278828573387727374*x(1)^3*x(3) - 0.00069288613969398227254714583978057*x(2)*x(3)^3 - 0.00052682852818064507260942264110781*x(1)^4*x(3) - 0.00051609977582134369811228680191562*x(2)^3*x(3) + 0.00021287226058633312852919061697321*x(2)^4*x(3) - 0.00013445172703319130391719227191061*x(1)^5 - 0.00024132867587522977004255153588019*x(2)^4 - 0.00020704747587133032205031213379698*x(3)^3 + 0.00011167455550850635681570111046312*x(3)^4 + 0.00040127472414558384983251926314551*x(1)*x(2)*x(3)^2 - 0.00021649604626866603140911138325464*x(1)*x(2)^2*x(3) + 0.00020608708784342066877570687211119*x(1)*x(2)^3*x(3) - 0.0011720311865741628309933730633929*x(1)^3*x(2)*x(3) - 0.00098316539033893590726620459463447*x(1)^2*x(2)^2*x(3); 
dx(2,1) = 0.0011317343059917828185234611737542*x(1)^2*x(2)^3 - 0.00024756957376509403090381056244951*x(1)^2*x(2)^2 + 0.00080955811953142831072227636468597*x(1)^2*x(3)^2 - 0.001348263315023023878325147961732*x(1)^3*x(2)^2 + 0.00033285791293075073227214488724712*x(1)^2*x(3)^3 + 0.00038900958090970494396287904237397*x(1)^3*x(3)^2 + 0.00016196806207777192376795483141905*x(2)^2*x(3)^2 + 0.00038832834796859172499239321041387*x(2)^2*x(3)^3 - 0.00042833367007794054259761651337612*x(2)^3*x(3)^2 - 0.00056748363721703665163431651308201*x(1)*x(2)^3 + 0.00011967817636705790906859192546108*x(1)*x(3)^2 - 0.0001901674501108285841866063492489*x(1)^3*x(2) + 0.00034168679695933956708131518098526*x(1)*x(2)^4 + 0.00078514990463873779447112610796466*x(1)*x(3)^3 + 0.00014313874311272511974379995081108*x(2)*x(3)^2 + 0.00023956671081959424185470197699033*x(1)^3*x(3) - 0.00082198902810193263945848229923286*x(1)^4*x(2) - 0.00013088910057648672768948472366901*x(1)*x(3)^4 + 0.0017131094754647691047466651070863*x(2)*x(3)^3 + 0.00031577732232668243028683718875982*x(1)^4*x(3) - 0.00013377229610345153787420713342726*x(2)^3*x(3) - 0.00023133388506677010632017754687695*x(2)*x(3)^4 + 0.00045787712063038998877573249046691*x(2)^4*x(3) + 0.0014643360876116506830157959484495*x(1)^5 - 0.0010759451939679198773092139163055*x(2)^4 + 0.0001911888501627545533523289122968*x(3)^3 - 0.00025307432821103370557125344930682*x(2)^5 + 0.00055721049697254887433928161044605*x(3)^4 + 0.00038667246550949663230767328059301*x(1)*x(2)*x(3)^2 + 0.00020896105815695897867101393785561*x(1)*x(2)*x(3)^3 - 0.00052137960131626304161045482032932*x(1)*x(2)^3*x(3) - 0.00060014009677067381431925241486169*x(1)^3*x(2)*x(3) + 0.00038870430498477714920113612606656*x(1)*x(2)^2*x(3)^2 - 0.00075323300072149823591871609096415*x(1)^2*x(2)*x(3)^2 - 0.00059113425555923360121823861845769*x(1)^2*x(2)^2*x(3); 
dx(3,1) = 0.00011139806290258458254527340614004*x(1)^3*x(3)^2 - 0.00048122246910592414437246588931885*x(1)^2*x(3)^2 - 0.00031519970027338306550745983258821*x(1)^2*x(2)^3 + 0.00060154606007845057291660850751214*x(2)^2*x(3)^2 + 0.00028318344455952049187885677383747*x(2)^3*x(3)^2 - 0.00017018299276375103978864444798091*x(1)*x(3)^2 - 0.00050391301062879811922812223201618*x(1)*x(2)^4 - 0.001984485603697017097601928981021*x(1)*x(3)^3 - 0.00041980956617648779172213835408911*x(1)^4*x(2) - 0.00056851180845751692061185167403892*x(2)*x(3)^3 + 0.00011198533399300880653015610732837*x(1)^4*x(3) - 0.00015728972778666916454426427662838*x(2)^3*x(3) + 0.00026971508813589117892206559190527*x(2)^4*x(3) + 0.00016739559778261581257652323984075*x(1)^5 - 0.00014151752785862559136376148671843*x(2)^4 - 0.00014804550860705867343369845912093*x(3)^3 + 0.0003538040161218947154964098444907*x(2)^5 - 0.00017219730205467187467149869917193*x(3)^4 + 0.00027637654889012264192160728271119*x(1)*x(2)*x(3)^2 - 0.00054530885143833085493270118604414*x(1)*x(2)^3*x(3) - 0.00073733816564069964982763849548064*x(1)^3*x(2)*x(3) + 0.00042692364538671201401598409574945*x(1)*x(2)^2*x(3)^2 + 0.00032832245469843757135208761610556*x(1)^2*x(2)*x(3)^2 - 0.00074235606241712659425502351950854*x(1)^2*x(2)^2*x(3); 
end