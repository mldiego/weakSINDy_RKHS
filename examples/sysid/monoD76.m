function dx = monoD(t,x) 
dx(1,1) = 0.12185325963038451391184935346246*x(1) + 0.25075540877219282265286892652512*x(2) + 0.058918733810614298818109091371298*x(3) + 0.04832287844634919338204781524837*x(1)^2*x(2)^2 + 0.0015412017010916212456095308880322*x(1)^2*x(2)^3 + 0.0061976637097211906279881077352911*x(1)^2*x(3)^2 + 0.0022618403194751834917042287997901*x(1)^3*x(2)^2 - 0.0001774908229437155426921890466474*x(1)^2*x(3)^3 + 0.001514919697301442269576909893658*x(1)^3*x(3)^2 + 0.02395529948845975809490482788533*x(2)^2*x(3)^2 - 0.0005495791468796973333610367262736*x(2)^2*x(3)^3 + 0.00029747688173398234567912368220277*x(2)^3*x(3)^2 + 0.15204799814728175988420844078064*x(1)*x(2) + 0.36367253691389578307280316948891*x(1)*x(3) + 0.83258685602902460232144221663475*x(2)*x(3) - 0.14468596291067115089390426874161*x(1)*x(2)^2 - 0.042904427806149669777369126677513*x(1)^2*x(2) - 0.14567793502484960299625527113676*x(1)*x(2)^3 - 0.234832265969600939570227637887*x(1)*x(3)^2 + 0.046147858690908094558835728093982*x(1)^2*x(3) + 0.12604831243325520517828408628702*x(1)^3*x(2) - 0.0018244339348227889274767221650109*x(1)*x(2)^4 + 0.01487124056479594003121746936813*x(1)*x(3)^3 + 0.031485870706553953368711518123746*x(2)*x(3)^2 - 0.0228498380464898787067795637995*x(1)^3*x(3) - 0.30150625460788660348043777048588*x(2)^2*x(3) - 0.00054801873429521918978934991173446*x(1)^4*x(2) - 0.00029605165654267295849422225728631*x(1)*x(3)^4 - 0.0036990054010397344086413795594126*x(2)*x(3)^3 + 0.00056229030517584988047019578516483*x(1)^4*x(3) - 0.016957790916404036352105322293937*x(2)^3*x(3) - 0.00051357290673925692914281171397306*x(2)^4*x(3) + 0.063411794519481645693304017186165*x(1)^2 - 0.071132028540674241412489209324121*x(1)^3 + 0.2568014870017805151292122900486*x(2)^2 - 0.022268909128534630781359737738967*x(1)^4 + 0.036694038961606167958962032571435*x(2)^3 + 0.073119177184551631398790050297976*x(3)^2 - 0.0014246941240407817730329043115489*x(1)^5 + 0.046823232015363203117885859683156*x(2)^4 - 0.013340636955309292943638865835965*x(3)^3 + 0.00021766568317510559538163761317264*x(2)^5 + 0.00073510853878833160024441895075142*x(3)^4 - 0.048317156953359585713769774883986*x(1)*x(2)*x(3)^2 + 0.096187773192554004708654247224331*x(1)*x(2)^2*x(3) - 0.073674758735108980545192025601864*x(1)^2*x(2)*x(3) + 0.0012502258664295329282367674750276*x(1)*x(2)*x(3)^3 + 0.0011331103769891548438408790389076*x(1)*x(2)^3*x(3) - 0.006775024658001171928844996728003*x(1)^3*x(2)*x(3) - 0.0021045744148859846234245196683332*x(1)*x(2)^2*x(3)^2 + 0.0010939012648842894037670703255571*x(1)^2*x(2)*x(3)^2 + 0.0032182989279365692425471934257075*x(1)^2*x(2)^2*x(3) + 0.44847124400547500044922344386578*x(1)*x(2)*x(3) + 0.009640135966831664404708135407418; 
dx(2,1) = 0.23298235111550980036554392427206*x(1) + 0.45990860899360086477827280759811*x(2) + 0.24891889022666191522148437798023*x(3) - 0.0073259058687149547495209844782948*x(1)^2*x(2)^2 + 0.0063991574762756187055856571532786*x(1)^2*x(3)^2 + 0.0021335723047255505946395715000108*x(1)^3*x(2)^2 - 0.0001744993506975212049781021050876*x(1)^2*x(3)^3 + 0.00038248581533084635708519272156991*x(1)^3*x(3)^2 + 0.018399766053729393888716003857553*x(2)^2*x(3)^2 - 0.00058740019122516251570687018102035*x(2)^2*x(3)^3 - 0.00059490618026303287990685930708423*x(2)^3*x(3)^2 + 0.33113373346606067570974119007587*x(1)*x(2) + 0.82828645545373547065537422895432*x(1)*x(3) + 1.5925099754324492096202448010445*x(2)*x(3) + 0.071381894927384337279363535344601*x(1)*x(2)^2 - 0.24022375155334430019138380885124*x(1)^2*x(2) + 0.0030777754252624411890337796648964*x(1)*x(2)^3 - 0.025865401006452515275668702088296*x(1)*x(3)^2 + 0.043380395458612497350259218364954*x(1)^2*x(3) - 0.105800126472146871492441277951*x(1)^3*x(2) + 0.00031098127900042182858442174619995*x(1)*x(2)^4 - 0.00035439329863468094217182624561246*x(1)*x(3)^3 - 0.15362363293357361726521048694849*x(2)*x(3)^2 - 0.00045212848251302562019304787099827*x(1)^3*x(3) - 0.1073466826468916224257554858923*x(2)^2*x(3) - 0.00011737609529482995185389881953597*x(1)^4*x(2) + 0.0048746608928631474100257037207484*x(2)*x(3)^3 + 0.00055295090800733248670439934358001*x(1)^4*x(3) - 0.0067391262006051277921869768761098*x(2)^3*x(3) - 0.00056869261471348409031634218990803*x(2)^4*x(3) + 0.15168933648914162404253147542477*x(1)^2 - 0.16205297600501467059075366705656*x(1)^3 + 0.77324248345564683404518291354179*x(2)^2 - 0.0049862824590194421148225956130773*x(1)^4 - 0.069206136783179772464791312813759*x(2)^3 + 0.09787985093893780685903038829565*x(3)^2 - 0.00044096020460471230251187080284581*x(1)^5 + 0.014413308513988809522743395064026*x(2)^4 - 0.020813109350456215906888246536255*x(3)^3 - 0.000625688303561799941121535084676*x(2)^5 + 0.0012803681128228383556688640965149*x(3)^4 - 0.0020488856734179883289925783174112*x(1)*x(2)*x(3)^2 + 0.032029294532598839850834337994456*x(1)*x(2)^2*x(3) + 0.04762182046201246521377470344305*x(1)^2*x(2)*x(3) + 0.00023981367366746053626513912604423*x(1)*x(2)*x(3)^3 + 0.00044342550380555012523586810857523*x(1)*x(2)^3*x(3) + 0.0022816063979123057947617780882865*x(1)^3*x(2)*x(3) - 0.0015966342800268762402993161231279*x(1)^2*x(2)*x(3)^2 + 0.00080779085028914732191651637549512*x(1)^2*x(2)^2*x(3) - 0.29590461402602841189946047961712*x(1)*x(2)*x(3) + 0.057886612870859721624583471566439; 
dx(3,1) = 0.00059760869039005015679322241339833*x(1)^2*x(3)^3 - 0.028188581396165091064176522195339*x(2) - 0.16939274422867356406641192734241*x(3) - 0.05787396657420629253465449437499*x(1)^2*x(2)^2 - 0.0013043210314052089415781665593386*x(1)^2*x(2)^3 - 0.031020064988133810857107164338231*x(1)^2*x(3)^2 - 0.0034819456762309464181726070819423*x(1)^3*x(2)^2 - 0.013847287258167639834027795586735*x(1) - 0.00032618568363812494581566170381848*x(1)^3*x(3)^2 - 0.041570308113684006912080803886056*x(2)^2*x(3)^2 + 0.0011799275889246008119926045765169*x(2)^2*x(3)^3 - 0.00047859018464452285357424443645868*x(2)^3*x(3)^2 + 0.018449302149157631447451421990991*x(1)*x(2) + 0.00041234371355702847949942224659026*x(1)*x(3) - 0.014833904414230048018907837104052*x(2)*x(3) - 0.040820448911802031943807378411293*x(1)*x(2)^2 - 0.060267647396123891212482703849673*x(1)^2*x(2) + 0.16865190552954345548641867935658*x(1)*x(2)^3 + 0.053961532275998536078986944630742*x(1)*x(3)^2 + 0.36634600938242556367185898125172*x(1)^2*x(3) - 0.15167744904212554502009879797697*x(1)^3*x(2) + 0.001753807241861249366365882451646*x(1)*x(2)^4 - 0.0061381801609421415832912316545844*x(1)*x(3)^3 - 0.0080727273611316263668413739651442*x(2)*x(3)^2 - 0.019994102233901855925068957731128*x(1)^3*x(3) + 0.31518545469066339137498289346695*x(2)^2*x(3) + 0.0039956136251424467786819150205702*x(1)^4*x(2) + 0.00017198624358835679082346814539051*x(1)*x(3)^4 + 0.0014265266781416929831038942211308*x(2)*x(3)^3 - 0.0033775575200083451932187017519027*x(1)^4*x(3) + 0.029042932820814826300193089991808*x(2)^3*x(3) + 0.0010425512793574842618227194179781*x(2)^4*x(3) + 0.066848119296068375660979654639959*x(1)^2 + 0.14829237021854169142898172140121*x(1)^3 + 0.16779090101897509157424792647362*x(2)^2 + 0.10157971477465821408259216696024*x(1)^4 + 0.023887276050409411709551932290196*x(2)^3 - 0.50910286075702515518059954047203*x(3)^2 + 0.00099005136449914488139256718568504*x(1)^5 - 0.056950472319584832803229801356792*x(2)^4 + 0.03828393307371413811779348179698*x(3)^3 - 0.00038680336345570109912728185008746*x(2)^5 - 0.0012844135051759408838734088931233*x(3)^4 + 0.072006097903241084168257657438517*x(1)*x(2)*x(3)^2 - 0.10279641069408285147801507264376*x(1)*x(2)^2*x(3) + 0.093624896839571647433331236243248*x(1)^2*x(2)*x(3) - 0.0017454950263939839061322345514782*x(1)*x(2)*x(3)^3 - 0.0011689855855576691823216606280766*x(1)*x(2)^3*x(3) + 0.0085805320921075178830506047233939*x(1)^3*x(2)*x(3) + 0.0021994339462465539725144481053576*x(1)*x(2)^2*x(3)^2 - 0.0018294376371803533487536697066389*x(1)^2*x(2)*x(3)^2 - 0.0055892438396734078764893638435751*x(1)^2*x(2)^2*x(3) - 0.57893871898988891189219430088997*x(1)*x(2)*x(3) - 0.031498992705380146617244463413954; 
end