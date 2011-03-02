{-# LANGUAGE TemplateHaskell, Rank2Types #-}
module Hyperbolic (
  hyperDist,
  arcLength,
  makeHyperArc,
  arcFrom,
  reflectPtThrough,
  reflectThrough
) where

import Data.Ratio ((%))
import Math (coshm1, squareEqSolutions, infinity, nan)
import Geometry (Point, Vector, fixAngle, (<.>), (</>), reflectVectorAgainst, dirVec)
import Primitives (Arc(..), arcNormals, arcFromPt, arcToPt)
import Data.SG (toAngle, mag, magSq, iso, fromPt, dotProduct, perpendicular2, distFrom, plusDir, minusDir, makeRel2, Point2'(..), unitVector, Pair(..), origin)
import Data.Ord (comparing)
import Data.List (find, minimumBy)
import Data.Maybe (fromJust)

import Test.QuickCheck
import Test.QuickCheck.All (forAllProperties)

import Debug.Trace


hyperDist :: Point -> Point -> Double
hyperDist u v | u2m1 > eps || v2m1 > eps = nan
              | q < eps   = infinity
              | otherwise = acosh (1 + delta u v)
  where
    eps = 1e-9
    u2m1 = magSq u - 1
    v2m1 = magSq v - 1
    q = u2m1*v2m1
    delta u v = 2 * (magSq $ u `fromPt` v) / q

arcLength :: Arc -> Double
arcLength arc = hyperDist (arcFromPt arc) (arcToPt arc)

makeHyperArc :: Point -> Point -> Arc
makeHyperArc u v  | abs q < eps = L { linePt1 = u, linePt2 = v }
                  | otherwise   = A { center = iso c, radius = r, fromA = a1fxd, toA = a2fxd }
  where
    eps = 1e-9
    q = (iso v) `dotProduct` (perpendicular2 (iso u)) :: Double
    c :: Point
    c = iso $ ((perpendicular2 (iso u)) <.> (magSq v + 1) - (perpendicular2 (iso v)) <.> (magSq u + 1)) </> (2*q)
    r = u `distFrom` c
    a1 = fixAngle $ toAngle $ u `fromPt` c
    a2 = fixAngle $ toAngle $ v `fromPt` c
    da = a2-a1
    (a1fxd, a2fxd)  | da < -pi  = (a1, a2 + 2*pi)
                    | da > pi   = (a1 + 2*pi, a2)
                    | otherwise = (a1, a2)

arcFrom :: Point -> Vector -> Maybe Arc
arcFrom u dir | abs uPerpDir < eps = Just $ L { linePt1 = u, linePt2 = lineEndPt }
              | otherwise = do ngonAngle <- return $ toAngle cosSinBeta
                               return $ A { center = c, radius = r, fromA = alpha, toA = ngonAngle }
  where
    eps = 1e-9
    -- v = u*k
    -- k^2 (1+bigU) - 2k + 1 - bigU / u2 == 0
    lineCoefSolutions = squareEqSolutions ((bigU+1), -2, 1-bigU/u2)
    lineEndPtSolutions = map (u <.>) lineCoefSolutions :: [Point]
    lineEndPt :: Point
    lineEndPt | null lineEndPtSolutions = u
              | otherwise = head $ filter (\v -> (v `fromPt` u) `dotProduct` dir > 0) lineEndPtSolutions
    
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
    alpha = toAngle $ iso $ u `fromPt` c :: Double

    -- various intermediate values used to compute the circle center/arc angles
    u2, bigD, bigU, bigR, bigX, bigY :: Double
    u2 = magSq u
    bigD = coshm1 (mag dir) / 2
    bigU = bigD * (1 - u2)
    bigR = let r2 = r^2 in bigU * (1 - magSq c - r2) / (2*r2) - 1
    Pair (bigX, bigY) = iso $ (c <.> (bigU/r)) `minusDir` (unitVector $ u `fromPt` c)

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

reflectPtThrough :: Point -> Arc -> Point
reflectPtThrough pt mirror = arcToPt imageArc
  where
    c = center mirror
    negMirrorDir = unitVector $ pt `fromPt` c
    mirrorPt = c `plusDir` (negMirrorDir <.> (radius mirror))
    mirrorNormal = perpendicular2 negMirrorDir
    preimageArc = makeHyperArc pt mirrorPt
    preimageArcNormal = negate $ snd $ snd $ arcNormals preimageArc
    reflectedNormal = preimageArcNormal `reflectVectorAgainst` mirrorNormal
    imageArc | arcLength preimageArc > 0 = fromJust $ arcFrom mirrorPt (reflectedNormal <.> (arcLength preimageArc)) -- FIXME: fromJust
             | otherwise = preimageArc

reflectThrough :: Arc -> Arc -> Arc
reflectThrough preimage mirror = makeHyperArc reflectedFrom reflectedTo
  where
    reflectedFrom = (arcFromPt preimage) `reflectPtThrough` mirror
    reflectedTo = (arcToPt preimage) `reflectPtThrough` mirror


----- Tests -----

arbFractional :: Fractional a => a -> a -> Gen a
arbFractional from to = sized $ \n -> do
    let maxN = ((toInteger n) * precision) `max` 1
    a <- choose (0, maxN)
    return $ fromRational (a % maxN) * (to-from)+from
  where precision = 9999999999999 :: Integer

arbAngle :: Floating a => Gen a
arbAngle = do
    m <- fmap fromIntegral $ choose (1, q)
    a <- fmap (*m) $ arbFractional 0 (twopi/(fromIntegral q))
    elements [-a, a]
  where
    q = 8*3*5 :: Int
    twopi = 2*pi

arbVector :: Double -> Double -> Gen Vector
arbVector minL maxL = do
    r <- arbFractional minL maxL
    a <- arbAngle
    return $ dirVec a <.> r

arbPoint :: Double -> Double -> Gen Point
arbPoint minL maxL = fmap (origin `plusDir`) $ arbVector minL maxL

arbHyperPoint, arbFiniteHyperPoint :: Gen Point
arbHyperPoint = arbPoint 0 1

arbFiniteHyperPoint = arbHyperPoint `suchThat` (\p -> magSq p < 1)
    
arbHyperArc, arbFiniteHyperArc :: Gen Arc
arbHyperArc = sized $ \n -> do
    c <- arbPoint 1 (fromIntegral $ n+1 `max` 10) `suchThat` (\p -> magSq p - 1 > eps)
    let r2 = magSq c - 1
    let r = sqrt r2
    -- now we need to find the range of valid angles of the arc
    let a = (r2 - 1 + magSq c) / (2 * mag c)
        alpha = acos $ a / r
        zeroDirAngle = toAngle $ origin `fromPt` c
        minAngle = zeroDirAngle - alpha
        maxAngle = zeroDirAngle + alpha
    fromAngle <- arbFractional minAngle maxAngle
    toAngle <- arbFractional minAngle maxAngle
    return $ A { center = c, radius = r, fromA = fromAngle, toA = toAngle }
  where
    eps = 1e-6

arbFiniteHyperArc = arbHyperArc `suchThat` (\a -> magSq (arcFromPt a) < 1 && magSq (arcToPt a) < 1)


prop_HyperDistInf = forAll (do
    u <- arbHyperPoint
    v <- arbPoint 1 1
    elements [(u,v),(v,u)]) $ \(u,v) -> isInfinite $ hyperDist u v 

prop_HyperDistNaN = forAll (sized $ \n -> do
    u <- arbHyperPoint
    v <- arbPoint 1 (fromIntegral n `max` 10) `suchThat` (\p -> magSq p - 1 > eps)
    elements [(u,v),(v,u)]) $ \(u,v) -> isNaN $ hyperDist u v 
  where
    eps = 1e-6

-- Failing:
-- - A {center = Point2 (-0.38701796244697073,46.772774248781474), radius = 46.763684562176365, fromA = -1.5595117191119663, toA = -1.5457892887965283}
--   Point2 (5.3359510955077634e-2,0.1893282277516213)
-- - A {center = Point2 (-7.026560251697408e-3,11.42896537535343), radius = 11.385134998039177, fromA = -1.6042656287329964, toA = -1.5008069125600594}
--   Point2 (2.6284348017659494e-3,-3.328187301160707e-3)
-- - A {center = Point2 (9.005954122824345,-1.674365647430688e-2), radius = 8.950278767303784, fromA = 3.1375978573775174, toA = 3.2429671227072405}
--   Point2 (0.26218718409077446,8.660415731694505e-3))
-- - A {center = Point2 (-1.5213160721601204,4.758137009633355), radius = 4.894310001813876, fromA = -1.3809081404953094, toA = -1.264400643335706}
--   Point4 (-7.157618553445397e-2,-0.838403729165129)
prop_ptDoubleReflect = forAll (do
    mirror <- arbHyperArc
    p <- arbFiniteHyperPoint
    return (mirror, p)) $ \(mirror, p) ->
      let reflect = (`reflectPtThrough` mirror)
        in hyperDist p (reflect $ reflect p) < eps
  where
    eps = 1e-6

allTests = $(forAllProperties) quickCheckResult
