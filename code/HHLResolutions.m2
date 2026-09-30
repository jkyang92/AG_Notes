newPackage(
    "HHLResolutions",
    Version => "0.3",
    Date => "July 16, 2026",
    Authors => {
        { Name => "Jay Yang"
        , Email => "jay.k.yang@vanderbilt.edu"}
        },
    Headline => "Code to work with HHL and related resolutions",
    AuxiliaryFiles => true,
    PackageExports => {"NormalToricVarieties","Complexes"},
    PackageImports => {"SimplicialComplexes"},
--    DebuggingMode => true,
    HomePage => "https://github.com/jkyang92/HHLResolutions/"
    )

export {
    "sliceByHyperplanes",
    "makeResolutionTable",
    "makeHHLPolytopes",
    "makeHHLPolytopesRelative",
    "makeHHLResolution",
    "hhlPolytopes",
    "hhlResolution",
    "hhlVectors",
    "toFacesByDimension",
    "toricSemigroupGens",
    "hhlModuleGens",
    "andersonDiagonalModuleGens",
    "andersonDiagonalModuleVertices",
    "andersonVertexToExponent",
    "andersonLaurentModule",
    "andersonModule",
    "hhlLaurentModule",
    "hhlModule",
    "hyperplaneStratificationPolytopes",
    "gensToLaurentModule",
    "gensToToricModule",
    "andersonDiagonalResolution",
    "bondalThomsenStrata",
    "makeResolution",
    "lineBundleBondalThomsenMonad",
    }

--a deprecation helper for all of the renaming
--provide symbols!
addDeprecatedName = (oldName,newName) -> (
    -- compute these values immediately so that errors trigger in installPackage/needsPackage
    -- instead of when the function is run
    f := value newName;
    errStr := "The function " | toString oldName | " has been renamed, please use " | toString newName | " instead";
    oldName <- (args -> (
            printerr errStr;
            f args
            )))


load "./HHLResolutions/resolution_tools.m2"
load "./HHLResolutions/hhl.m2"
load "./HHLResolutions/anderson.m2"
load "./HHLResolutions/modules.m2"
load "./HHLResolutions/monads.m2"

load "./HHLResolutions/tests.m2"
load "./HHLResolutions/doc.m2"

end--

restart
installPackage("HHLResolutions", RemakeAllDocumentation=>  true)
installPackage"HHLResolutions"
needsPackage "HHLResolutions"
check HHLResolutions
viewHelp HHLResolutions

X = weightedProjectiveSpace {2,6,15,30}

gcd {6,10,30,45}

X = weightedProjectiveSpace {42,195,286,385}

factor 42
factor 195
X = weightedProjectiveSpace {30,154,429,455}

{30,154,429,455} / factor

p1 = 2
p2 = 3
p3 = 5
p4 = 7
p5 = 11
p6 = 13

gcd {p1*p2*p3,p1*p4*p5,p2*p4*p6,p3*p5*p6}
X = weightedProjectiveSpace {2*3*5,2*3*7,2*5*7,3*5*7}

Y = normalToricVariety({{0}},{{0}})
phi = map(X,Y,0);
    S = ring X;
    C = makeHHLResolution(X,matrix phi)
sort (degrees C.dd_1)_0

matrix phi

C1 = makeHHLResolution(X,matrix phi)
C2 = makeHHLResolution(X,matrix {{},{},{}})

(cells,verts,g,box) = makeHHLPolytopes(X,matrix {{},{},{}})

vertices convexHull cells

pruneComplex C2

degrees prune HH_0 C2
gcd {42,195,286,385}


sort degrees C_0
sort degrees C_1

sort (degrees C1.dd_0)_1
sort (degrees C2.dd_1)_0

dim Y

X = toricProjectiveSpace 3
Y = X ** X;
phi = diagonalToricMap X
S = ring Y;

hhlModuleGens(phi)
andersonDiagonalModuleGens(X)
andersonDiagonalModuleVertices(X)

X = hirzebruchSurface 2
Y = X ** X;
phi = diagonalToricMap X
S = ring Y;
G := hhlModuleGens(phi)
G' := andersonDiagonalModuleGens(X)
V := andersonDiagonalModuleVertices(X)
apply(V,v -> andersonVertexToExponent(X,v))

netList pairs bondalThomsenStrata rays X

ML := module G_0
ML_7
ML_6

X = normalToricVariety ({{1,0},{1,1},{0,1},{-1,1},{0,-1}},{{0,1},{1,2},{2,3},{3,4},{4,0}})
Y = X ** X;
phi = diagonalToricMap X
S = ring Y;
G := hhlModuleGens(phi)
G' := andersonDiagonalModuleGens(X)
ML := module G'_0

(ML_1+ML_9-(ML_4+ML_6))
G'_1

ML_0+ML_1+ML_8-(ML_5+ML_6+ML_3)

C := makeHHLResolution(Y,matrix phi);

M2 := trim HH_0 C

prune M1
prune M2
prune comodule M1
prune M2


X = normalToricVariety ({{1,0},{1,1},{0,1},{-1,1},{0,-1}},{{0,1},{1,2},{2,3},{3,4},{4,0}})
Y = X**X
S = ring Y

matrix rays Y
degreeTable := degrees S


degree (S_1*S_9) - degree (S_4*S_6)

C := andersonDiagonalResolution X

C.dd^2 == 0

degrees C_0
degrees C_1



primaryDecomposition ann HH_0 andersonDiagonalResolution X

decompose ideal Y
