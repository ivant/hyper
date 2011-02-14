module Hyperbolic where

import Graphics.Rendering.Cairo
import Control.Arrow hiding ((<+>))
import Control.Applicative
import Data.Ord
import Data.List (find, minimumBy)
import Data.SG hiding (lineTo)
import Data.SG.Geometry hiding (lineTo)
import Data.SG.Geometry.TwoDim
import Data.SG.Vector
import Data.SG.Vector.Basic
import Text.Printf
import Data.Maybe (fromJust)

import Debug.Trace

type Point = Point2' Double
data Arc = A { center :: Point, radius :: Double, fromA :: Double, toA :: Double } deriving Show

type Vector = Rel2' Double

(<+>), (<->) :: VectorNum v => v Double -> v Double -> v Double
a <+> b = fmapNum2 (+) a b
a <-> b = fmapNum2 (-) a b

(<.>), (</>) :: VectorNum v => v Double -> Double -> v Double
v <.> s = (*s) `fmapNum1` v
v </> s = (/s) `fmapNum1` v

infixl 6 <+>, <->
infixl 7 <.>, </>

-- expm1 x = e^x - 1, accurate for small x
expm1 x | u == 1    = x
        | um1 == -1 = -1
        | otherwise = um1*x/(log u)
  where
    u = exp x
    um1 = u - 1
-- coshm1 x = cosh(x) - 1, accurate for small x
-- cosh(x) = (e^x + e^(-x))/2-1 =
--         = (e^x - 2 + e^(-x))/2 =
--         = (e^(2x) - 2e^x + 1)/(2e^x) =
--         = (e^x - 1)^2 / (2e^x) =
--         = expm1(x)^2 / (2*(expm1(x)+1))
coshm1 :: Double -> Double
coshm1 x = let exm1 = expm1 x in exm1^2 / (2*(exm1+1))

dirVec :: Double -> Vector
dirVec a = makeRel2 (cos a, sin a)

rotateVec :: Double -> Vector -> Vector
rotateVec a v = makeRel2 (t `dotProduct` v, tp `dotProduct` v)
  where
    t = makeRel2 (cos a, -sin a)
    tp = makeRel2 (sin a, cos a)

hyperDist u v = acosh (1 + delta u v)
  where
    delta u v = 2 * (magSq $ u <-> v) / ((1 - magSq u)*(1 - magSq v))

arcLength arc = hyperDist (fst $ fst normals) (fst $ snd normals)
  where
    normals = arcNormals arc

arcFromPt, arcToPt :: Arc -> Point
arcFromPt arc = (center arc) `plusDir` ((dirVec $ fromA arc) <.> (radius arc))
arcToPt   arc = (center arc) `plusDir` ((dirVec $ toA arc) <.> (radius arc))

segmentToArc :: Point -> Point -> Arc
segmentToArc u v = A { center = iso c, radius = r, fromA = a1fxd, toA = a2fxd }
  where
    q = (iso v) `dotProduct` (perpendicular2 (iso u))
    c :: Point
    c = iso $ ((perpendicular2 (iso u)) <.> (magSq v + 1) - (perpendicular2 (iso v)) <.> (magSq u + 1)) </> (2*q)
    r = u `distFrom` c
    a1 = fixAngle $ toAngle $ u `fromPt` c
    a2 = fixAngle $ toAngle $ v `fromPt` c
    da = a2-a1
    (a1fxd, a2fxd)  | da < -pi  = (a1, a2 + 2*pi)
                    | da > pi   = (a1 + 2*pi, a2)
                    | otherwise = (a1, a2)
    {- RSI: remove
    (a1fxd, a2fxd) =
      if da < 0 then
          if da < -pi then
              (a1, a2+2*pi)
            else
              (a2, a1)
        else
          if da > pi then
              (a2, a1+2*pi)
            else
              (a1, a2)
    -}

fixAngle a  | a < -pi   = fixAngle (a+2*pi)
            | a > pi    = fixAngle (a-2*pi)
            | otherwise = a

arcNormals :: Arc -> ((Point, Vector), (Point, Vector))
arcNormals arc = ((p1,v1), (p2, v2))
  where
    p1 = arcFromPt arc
    p2 = arcToPt arc
    v1 = (perpendicular2 $ dirVec $ fromA arc) <.> signum a
    v2 = (perpendicular2 $ dirVec $ toA arc) <.> signum a
    a = fixAngle (toA arc - fromA arc)

arcFrom :: Point -> Vector -> Maybe Arc
arcFrom u dir | abs uPerpDir < eps = Nothing
              | otherwise = do ngonAngle <- return $ toAngle cosSinBeta
                               return $ A { center = c, radius = r, fromA = alpha, toA = ngonAngle }
  where
    eps = 1e-9
    -- center of the arc circle is defined by the following equations:
    -- 2 (u.perpDir) x0 = u1 (u.perpDir) + u2 (u.dir) - d2 =       u . ( u.perpDir, u.dir ) - d2
    -- 2 (u.perpDir) y0 = u2 (u.perpDir) - u1 (u.dir) + d1 = - perpU . ( u.perpDir, u.dir ) + d1
    c :: Point
    c = (Point2 (iso u `dotProduct` uDirV, negate $ perpU `dotProduct` uDirV) `plusDir` perpDir) </> (2*uPerpDir)
      where
        uDirV = makeRel2 $ (uPerpDir, iso u `dotProduct` dir) :: Vector
        perpU = perpendicular2 $ iso u                        :: Vector
    perpDir = perpendicular2 dir                              :: Vector
    uPerpDir = iso u `dotProduct` perpDir                     :: Double
    -- radius of the arc circle
    r = u `distFrom` c :: Double
    alpha = toAngle $ iso $ u <-> c :: Double

    -- various intermediate values used to compute the circle center/arc angles
    bigD, bigU, bigR, bigX, bigY :: Double
    bigD = coshm1 (mag dir) / 2
    bigU = bigD * (1 - magSq u)
    bigR = let r2 = r^2 in bigU * (1 - magSq c - r2) / (2*r2) - 1
    Pair (bigX, bigY) = iso $ c <.> (bigU/r) <-> (unitVector $ u <-> c)

    bigXSq = bigX^2
    bigYSq = bigY^2
    bigRSq = bigR^2

    -- sin and cos of the angle of the end-point of the arc, with the circle center as a reference point
    -- may contain illegal values (outside of [-1,1] range)
    cosBetaSolutions, sinBetaSolutions :: [Double]
    (cosBetaSolutions, sinBetaSolutions)  | bigXSq > bigYSq = (c2, s2)
                                          | otherwise       = (c1, s1)
      where
        s1 = squareEqSolutions (bigXSq+bigYSq, -2*bigR*bigY, bigRSq-bigXSq)
        c2 = squareEqSolutions (bigXSq+bigYSq, -2*bigR*bigX, bigRSq-bigYSq)

        c1  | bigXSq > 0  = map (\sinBeta -> (bigR - bigY*sinBeta)/bigX) s1
            | otherwise   = concatMap (\sinBeta -> let v = sqrt (1 - sinBeta^2) in [v, -v]) s1
        s2  | bigYSq > 0  = map (\cosBeta -> (bigR - bigX*cosBeta)/bigY) c2
            | otherwise   = concatMap (\cosBeta -> let v = sqrt (1 - cosBeta^2) in [v, -v]) c2

    -- same as {sin,cos}BetaSolutions, but the illegal pairs of values have been filtered out
    cosSinBetaValid :: [Vector]
    cosSinBetaValid = filter validCosSinPair $ map makeRel2 $ zip cosBetaSolutions sinBetaSolutions
      where
        validCosSinPair :: Vector -> Bool
        validCosSinPair ngonAngleV = abs (magSq ngonAngleV - 1) < eps && magSq (c `plusDir` (ngonAngleV <.> r)) <= 1


    -- end-point direction (from the arc circle center)
    cosSinBetas :: [Vector]
    cosSinBetas = filter (\ngonAngleV -> ((vdir ngonAngleV) `dotProduct` dir) >= 0) cosSinBetaValid
      where
        vdir ngonAngleV = (c `plusDir` (ngonAngleV <.> r)) `fromPt` u

    vs :: [Point]
    vs = map ((c `plusDir`) . (<.> r)) cosSinBetas

    cosSinBeta :: Vector
    cosSinBeta = snd $ minimumBy (comparing ((hyperDist u).fst)) $ zip vs cosSinBetas


    tr x = trace (show x) x
    -- solve a x^2 + b x + c == 0
    squareEqSolutions :: (Double, Double, Double) -> [Double]
    squareEqSolutions (a,b,c) | a == 0 = if b /= 0 then [-c / (2*b)] else []
                              | d > 0  = [parabolaTop + intersectionOffset, parabolaTop - intersectionOffset]
                              | d == 0 = [parabolaTop]
                              | d < 0  = []
                              | otherwise = error $ "d is NaN in a square equation"
      where
        d = b^2 - 4*a*c
        intersectionOffset = sqrt d / (2*a)
        parabolaTop = -b / (2*a)

setupPNG :: Int -> Render ()
setupPNG sz = do
    let sc = 1.3
    scale (fromIntegral sz / (2*sc)) (-fromIntegral sz / (2*sc))
    translate (sc) (-sc)
    setLineWidth $ ((2*sc) / fromIntegral sz)

reflectPtThrough :: Point -> Arc -> Point
reflectPtThrough pt mirror = arcToPt imageArc
  where
    (mirrorPt, mirrorNormal) = fst $ arcNormals mirror
    preimageArc = segmentToArc pt mirrorPt
    preimageArcNormal = negate $ snd $ snd $ arcNormals preimageArc
    reflectedNormal = preimageArcNormal `reflectVectorAgainst` mirrorNormal
    imageArc | arcLength preimageArc > 0 = fromJust $ arcFrom mirrorPt (reflectedNormal <.> (arcLength preimageArc)) -- FIXME: fromJust
             | otherwise = preimageArc

reflectThrough :: Arc -> Arc -> Arc
reflectThrough preimage mirror = segmentToArc reflectedFrom reflectedTo
  where
    reflectedFrom = (arcFromPt preimage) `reflectPtThrough` mirror
    reflectedTo = (arcToPt preimage) `reflectPtThrough` mirror

reflectVectorAgainst :: Vector -> Vector -> Vector
reflectVectorAgainst a b = c <.> 2 - au
  where
    au = unitVector a
    bu = unitVector b
    c = b <.> (au `dotProduct` bu)

drawHyper :: Render ()
drawHyper = do
    setSourceRGB 1 1 1
    paint

    setSourceRGB 0 0 0
    setDash [0.1, 0.1] 0.15
    moveTo 0 (-1.5)
    lineTo 0 1.5
    moveTo (-1.5) 0
    lineTo 1.5 0
    stroke
    arc 0 0 1 0 (2*pi)
    stroke

    setDash [] 0

    mapM_ (\a -> drawArc (0,0,1) a) ngon
    {-
    let bs = [ a `reflectThrough` m | a <- ngon, m <- ngon ]
    mapM (\b -> drawArc (0,0.5,0.7) b) bs
    -}
    mapM_ (\a -> drawArc (0,0.5,0.7) a) $ map (`reflectThrough` mirror) [arcToMirror]

    drawPoint (1,0,0) pt
    drawArc (0.7,0,0.5) mirror
--    drawArc (0.7,0.5,0) preimageArc
--    drawArc (0.7,0.5,0) imageArc

    return ()

  where
    arcToMirror = ngon !! 0
    mirror = ngon !! 1
    pt = arcToPt arcToMirror

    (mirrorPt, mirrorNormal) = fst $ arcNormals mirror
    preimageArc = segmentToArc pt mirrorPt
    preimageArcNormal = negate $ snd $ snd $ arcNormals preimageArc
    reflectedNormal = preimageArcNormal `reflectVectorAgainst` mirrorNormal
    imageArc | arcLength preimageArc > 0 = fromJust $ arcFrom mirrorPt (reflectedNormal <.> (arcLength preimageArc)) -- FIXME: fromJust
             | otherwise = preimageArc

    drawStraight (r,g,b) pt1 pt2 =
      let Pair (x0,y0) = iso pt1
          Pair (x1,y1) = iso pt2
        in do
          save
          setSourceRGB r g b
          moveTo x0 y0
          lineTo x1 y1
          stroke
          restore

    drawPoint (r,g,b) pt = do
      save
      let Pair (x,y) = iso pt
      setSourceRGB 1 0 0
      moveTo x y
      arc x y 0.01 0 (2*pi)
      stroke
      restore

    drawNormals a = do
      save
      setSourceRGB 1 0 0
      setDash [] 0
      let arrlen = 0.07
      let n = fst $ arcNormals a
      let Pair (x0,y0) = iso $ fst n
      let Pair (x1,y1) = iso $ (fst n) `plusDir` ((snd n) <.> arrlen)
      moveTo x0 y0
      lineTo x1 y1
      let n = snd $ arcNormals a
      let Pair (x0,y0) = iso $ fst n
      let Pair (x1,y1) = iso $ (fst n) `plusDir` ((snd n) <.> arrlen)
      moveTo x0 y0
      lineTo x1 y1
      stroke
      restore

    drawArcConstruction a = do
      save
      setLineWidth (1/(8*72))   -- FIXME calculate line width based on current scale
      setDash [0.01, 0.01] 0.015
      setSourceRGB 0.5 0.5 0.5
      arc (getX (center a)) (getY (center a)) (radius a) 0 (2*pi)
      stroke
      restore

    drawArc (r,g,b) a = do
      --drawArcConstruction a

      save
      setDash [] 0
      setSourceRGB r g b
      let arcF = if (fixAngle $ toA a - fromA a) > 0 then arc else arcNegative
      arcF (getX (center a)) (getY (center a)) (radius a) (fromA a) (toA a)
      stroke
      restore

      --drawNormals a

arcAngles :: Arc -> Arc -> Double
arcAngles a1 a2 = 180 * (acos $ (unitVector a2n) `dotProduct` (unitVector a1n)) / pi
  where
    a1n = snd $ snd $ arcNormals a1
    a2n = snd $ fst $ arcNormals a2

nextArc alpha len arc = arcFrom (fst n) dir
  where
    n = snd $ arcNormals arc
    dir = (rotateVec alpha (unitVector $ snd n)) <.> len

ngonAngle = pi/2
ngonSides = 7 :: Int
beta = 2*pi/(fromIntegral ngonSides)
len = acosh (1+2*(cos beta))

ngon = take ngonSides $ iterate (fromJust . nextArc ngonAngle len) firstArc
  where
    firstArc = fromJust $ arcFrom startP (dirVec (ngonAngle/2+pi/2) <.> len)
    startP = Point2 (r, 0)
    r = let t = tan (beta/2) in sqrt ((1-t)/(1+t))

{-
Just a1 = arcFrom startP (dirVec (ngonAngle/2+pi/2) <.> len)
Just a2 = nextArc ngonAngle len a1
Just a3 = nextArc ngonAngle len a2
Just a4 = nextArc ngonAngle len a3
Just a5 = nextArc ngonAngle len a4
as = [a1, a2, a3, a4, a5]
-}


main = do
  let sz = 600
  {-
  mapM_ (printf "%f\n" . (subtract $ len) . arcLength) as
  printf "a12: %f\n" $ arcAngles a1 a2
  printf "a23: %f\n" $ arcAngles a2 a3
  printf "a34: %f\n" $ arcAngles a3 a4
  printf "a45: %f\n" $ arcAngles a4 a5
  printf "d(a5,startP) = %f\n" $ (fst$snd$arcNormals a5) `hyperDist` startP
  -}
  withImageSurface FormatRGB24 sz sz $ \s -> do
    renderWith s $ (setupPNG sz >> drawHyper)
    surfaceWriteToPNG s "hyper.png"

  return ()
