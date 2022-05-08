function dx = monoD(t,x) 
dx(1,1) = 0.0001488582744489097553497458648053*x(2) - 0.0041905666321362389226123923435807*x(1)^2*x(2)^2 + 0.0049409511059765875984339800197631*x(1)^2*x(2)^3 - 0.0055007586745219327895028982311487*x(1)^2*x(3)^2 - 0.0026947683337357197785877360729501*x(1)^3*x(2)^2 + 0.00040715424965925839373426242673304*x(1)^3*x(3)^2 + 0.01764082791834553631815651897341*x(2)^2*x(3)^2 - 0.00063552353749518619707714606192894*x(2)^2*x(3)^3 - 0.00079740982904896728911126047023572*x(2)^3*x(3)^2 + 0.00059849568968828315007613127818331*x(1)*x(3) + 0.0013637113729825234287318380665965*x(2)*x(3) - 0.0025268692276751991698802157770842*x(1)*x(2)^2 - 0.0012412858513923819003821336082183*x(1)^2*x(2) - 0.00088441786086168594493983619031496*x(1)*x(2)^3 + 0.0036887828904395014717465528519824*x(1)*x(3)^2 - 0.00030224274753459789977227956114803*x(1)^2*x(3) - 0.0014327811009844371170629528933205*x(1)^3*x(2) + 0.00079589956242220871729386999504641*x(1)*x(2)^4 + 0.00093884104032759818636577620054595*x(1)*x(3)^3 + 0.0095292481844460041884303791448474*x(2)*x(3)^2 - 0.0048337993218519059723803366068751*x(1)^3*x(3) + 0.004663624080358808043911267304793*x(2)^2*x(3) - 0.0055529512198075892115411988925189*x(1)^4*x(2) - 0.00020334307218133895744927031046245*x(2)*x(3)^3 + 0.0010324720523648966974405993823893*x(1)^4*x(3) + 0.0027208912243854577184265508549288*x(2)^3*x(3) - 0.00043324989461773766308283484249841*x(2)^4*x(3) - 0.00068193843283037569591442661476322*x(1)^3 + 0.00011007462877606066786739802410011*x(1)^4 - 0.0061745118292071410337484849151224*x(2)^3 - 0.00066250160129699420252791242091917*x(3)^2 - 0.002704174579752915974495408590883*x(1)^5 + 0.0020219302049668819165617605904117*x(2)^4 - 0.0099338963176958117173853679560125*x(3)^3 - 0.00087152688692193613917424954706803*x(2)^5 + 0.00060845791453700837081441932241432*x(3)^4 - 0.0071199866091617991514794994145632*x(1)*x(2)*x(3)^2 - 0.0075378593617179134867001266684383*x(1)*x(2)^2*x(3) - 0.0085565859640084340753674041479826*x(1)^2*x(2)*x(3) - 0.00016248505168126037290221574949101*x(1)*x(2)*x(3)^3 - 0.0010949512296456465065830343519337*x(1)*x(2)^3*x(3) + 0.0013675644213064774845634019584395*x(1)^3*x(2)*x(3) + 0.001013266533428813787054423301015*x(1)*x(2)^2*x(3)^2 + 0.00032723947678919840242883765313309*x(1)^2*x(2)*x(3)^2 + 0.0032057839328993864569383731577545*x(1)^2*x(2)^2*x(3) + 0.0012011126584454689947278893669136*x(1)*x(2)*x(3); 
dx(2,1) = 0.00014774330469363716744624070997816*x(1) + 0.00028118409674926514441040126257576*x(2) - 0.0056044661000855811039400578010827*x(1)^2*x(2)^2 + 0.0066727109534854278649618208874017*x(1)^2*x(2)^3 + 0.013341082387711367118754424154758*x(1)^2*x(3)^2 - 0.00432905690416074406812185770832*x(1)^3*x(2)^2 - 0.00032562617631831658471242008090485*x(1)^2*x(3)^3 + 0.00028592280840938899544312334910501*x(1)^3*x(3)^2 + 0.022127477238825576932867988944054*x(2)^2*x(3)^2 - 0.00093889017460102586198900098679587*x(2)^2*x(3)^3 - 0.0011232276585815270664170384407043*x(2)^3*x(3)^2 + 0.0011857806145247895557304218527861*x(1)*x(3) + 0.0022819517136754008390653325477615*x(2)*x(3) - 0.0072842042084388936018513049930334*x(1)*x(2)^2 - 0.0037074535141767839263593486975878*x(1)^2*x(2) + 0.00036715460094199903195999468152877*x(1)*x(2)^3 + 0.0072450877098155785915878368541598*x(1)*x(3)^2 + 0.0017928364478803082704416738124564*x(1)^2*x(3) - 0.0022143484328966955843043251661584*x(1)^3*x(2) + 0.0011413713546852832791955734137446*x(1)*x(2)^4 + 0.0032010691013470982113631180254743*x(1)*x(3)^3 + 0.015902803351359295902511803433299*x(2)*x(3)^2 - 0.0040714190159105712041309743653983*x(1)^3*x(3) + 0.008475962250544810672181483823806*x(2)^2*x(3) - 0.0076699173822865063243625627364963*x(1)^4*x(2) - 0.00012964773169737986435734455881175*x(1)*x(3)^4 - 0.0011203834704196502514150779461488*x(2)*x(3)^3 + 0.0025698455666001152053468103986233*x(1)^4*x(3) + 0.0089862557988240610029606614261866*x(2)^3*x(3) - 0.00054186751086160445112227534991689*x(2)^4*x(3) - 0.001431777990784199872109638818074*x(1)^3 + 0.0018481077906726817161597864469513*x(1)^4 - 0.016386620001345164610029314644635*x(2)^3 - 0.00072170850520214546719444115296938*x(3)^2 - 0.00024013932791064340754871864191955*x(1)^5 + 0.0013119446735914142010415162076242*x(2)^4 - 0.01172660859237062425108888419345*x(3)^3 - 0.0012220172066776946451227559009567*x(2)^5 + 0.0010992616527860654684900509892032*x(3)^4 - 0.0098655495596204900721204467117786*x(1)*x(2)*x(3)^2 - 0.016109905560835358073745737783611*x(1)*x(2)^2*x(3) - 0.019287380382991159422090277075768*x(1)^2*x(2)*x(3) + 0.00010404385826864870345787039696006*x(1)*x(2)*x(3)^3 - 0.0019347912866272132248468551551923*x(1)*x(2)^3*x(3) + 0.0025579577363927796795906033366919*x(1)^3*x(2)*x(3) + 0.0015516097425456454317327370517887*x(1)*x(2)^2*x(3)^2 - 0.00060797371293641599976353973033838*x(1)^2*x(2)*x(3)^2 + 0.0045082616945206765990405983757228*x(1)^2*x(2)^2*x(3) + 0.0028964406557077104764630348654464*x(1)*x(2)*x(3); 
dx(3,1) = 0.0055290815892483635707321809604764*x(1)^2*x(3)^2 - 0.0072084673003010024672221334185451*x(1)^2*x(2)^3 - 0.0039661705161062599245269666425884*x(1)^2*x(2)^2 + 0.0063551715930278973587519431021065*x(1)^3*x(2)^2 - 0.00015499642882418895872831399174174*x(1)^2*x(3)^3 - 0.00036195203058897451597886174567975*x(1)^3*x(3)^2 + 0.0027040630352939132308165426366031*x(2)^2*x(3)^2 - 0.00058208327993236608222105132881552*x(2)^3*x(3)^2 + 0.00032747750436845635135796328540891*x(1)*x(2) + 0.00020670059421262765297910846129525*x(1)*x(3) + 0.0003848340872754651975640172167914*x(2)*x(3) + 0.00051241027120851878606799800763838*x(1)*x(2)^2 + 0.00017103553660685610893210650829133*x(1)^2*x(2) - 0.0051235917778837958280746534001082*x(1)*x(2)^3 + 0.00042619454628123332184941318701021*x(1)*x(3)^2 + 0.0011060160917706340910626749973744*x(1)^2*x(3) - 0.0015195680992590165914180033723824*x(1)^3*x(2) + 0.0039664381538440096619524410925806*x(1)*x(2)^4 - 0.00012955998374367116987571080244379*x(1)*x(3)^3 + 0.0005919669717334707925715520104859*x(2)*x(3)^2 - 0.0011769351096400093581451073987409*x(1)^3*x(3) + 0.0055257205175127666052503627724946*x(2)^2*x(3) - 0.0010358226219673305479318514699116*x(1)^4*x(2) - 0.006091130030526414884661789983511*x(2)*x(3)^3 + 0.0050904686339716676002353779040277*x(2)^3*x(3) + 0.00027705229936719755201579573622439*x(2)*x(3)^4 - 0.00053603539137014877269393764436245*x(2)^4*x(3) + 0.00014938387912000195001382962800562*x(1)^2 + 0.00069272810748877144959578799898736*x(2)^2 - 0.00083232091078366821079725923482329*x(1)^4 + 0.0014510713325837087950276327319443*x(2)^3 - 0.00084823517928112046604383067460731*x(3)^2 + 0.0044937933933439921929675620049238*x(1)^5 + 0.0028181571074528122267111029941589*x(2)^4 - 0.0074275410923068108104416751302779*x(3)^3 - 0.0010126681535878034878805920016021*x(2)^5 + 0.0069069035717550519848373369313776*x(1)*x(2)*x(3)^2 - 0.0011848700583858562396244451520033*x(1)*x(2)^2*x(3) + 0.00017541221919631366787939441564959*x(1)^2*x(2)*x(3) - 0.00054784001167162976742019964149222*x(1)*x(2)*x(3)^3 + 0.0024733647787318346900065080262721*x(1)*x(2)^3*x(3) - 0.0029372442084367733627914276439697*x(1)^3*x(2)*x(3) + 0.00092032940231479010151360853342339*x(1)*x(2)^2*x(3)^2 + 0.00039197029606552469260805082740262*x(1)^2*x(2)*x(3)^2 - 0.0023908209222689436046493938192725*x(1)^2*x(2)^2*x(3) + 0.0024739656644507945770783408079296*x(1)*x(2)*x(3); 
end