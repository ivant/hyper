module Hyperbolic (
  hyperDist,
  arcLength,
  makeHyperArc,
  arcFrom,
  reflectPtThrough,
  reflectThrough
) where

import Math (coshm1, squareEqSolutions)
import Geometry (Point, Vector, fixAngle, (<.>), (</>), reflectVectorAgainst)
import Primitives (Arc(..), arcNormals, arcFromPt, arcToPt)
import Data.SG (toAngle, mag, magSq, iso, fromPt, dotProduct, perpendicular2, distFrom, plusDir, minusDir, makeRel2, Point2'(..), unitVector, Pair(..))
import Data.Ord (comparing)
import Data.List (find, minimumBy)
import Data.Maybe (fromJust)

import Debug.Trace


hyperDist :: Point -> Point -> Double
hyperDist u v = acosh (1 + delta u v)
  where
    delta u v = 2 * (magSq $ u `fromPt` v) / ((1 - magSq u)*(1 - magSq v))

arcLength :: Arc -> Double
arcLength arc = hyperDist (fst $ fst normals) (fst $ snd normals)
  where
    normals = arcNormals arc

makeHyperArc :: Point -> Point -> Arc
makeHyperArc u v  | q == 0 = L { linePt1 = u, linePt2 = v }
                  | otherwise = A { center = iso c, radius = r, fromA = a1fxd, toA = a2fxd }
  where
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
    lineCoefSolutions = squareEqSolutions ((bigU+1)*(magSq dir), -2*(mag u)*(mag dir), u2-bigU)
    lineEndPtSolutions = map (\k -> dir <.> k) $ filter (>=0) lineCoefSolutions
    lineEndPt :: Point
    lineEndPt | null lineEndPtSolutions = u
              | otherwise = iso $ head lineEndPtSolutions
    
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
    (mirrorPt, mirrorNormal) = fst $ arcNormals mirror
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

