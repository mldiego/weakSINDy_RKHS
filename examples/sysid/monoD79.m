function dx = monoD(t,x) 
dx(1,1) = 9.9990665334971708944067358970642*x(2) - 9.9981456288605841109529137611389*x(1) + 0.00032623387652885860177320864750072*x(3) - 0.00025035411230661663850582954182755*x(1)*x(3) + 0.00012075782679195345981071341157076*x(2)*x(3) - 0.0011615428965692231599859951529652; 
dx(2,1) = 27.997113805715343914926052093506*x(1) - 0.99863493543671211227774620056152*x(2) - 0.00016745793121292207317196698568296*x(3) - 0.99967099172829421149799600243568*x(1)*x(3) - 0.00014550467562862712256332997640129*x(2)*x(3) + 0.00011182372264932627370370710195857*x(1)^2 + 0.00061393479579763265974179375916719; 
dx(3,1) = 0.0015226107649470854710216372041032*x(1) - 0.00064030553938687617687719466630369*x(2) - 2.6665228603719697275664657354355*x(3) + 1.0000023017380499368300661444664*x(1)*x(2) - 0.00019539428437620465217605669749901*x(1)*x(3) - 0.00072227701529203880426166506367736; 
end