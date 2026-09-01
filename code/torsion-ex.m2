needsPackage "NormalToricVarieties"

ray = {{2,1},{0,1},{-2,-1},{0,-1}}
con = {{0,1},{1,2},{2,3},{3,0}}
X = normalToricVariety(ray,con)

isWellDefined X
isSmooth X
isSimplicial X

classGroup X
fromPicToCl X

S = ring X
deg = transpose matrix degrees S
