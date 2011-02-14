module Geometry (
  Point,
  Vector,
  (<.>), (</>),
  dirVec,
  fixAngle,
  rotateVec,
  reflectVectorAgainst
) where

import Data.SG (VectorNum(..), dotProduct, Point2', Rel2', makeRel2, unitVector)

type Point = Point2' Double
type Vector = Rel2' Double

-- | Scale vector by the given value
(<.>) :: VectorNum v => v Double -> Double -> v Double
v <.> s = (*s) `fmapNum1` v

-- | Scale vector by the reciprocal of the given value
(</>) :: VectorNum v => v Double -> Double -> v Double
v </> s = (/s) `fmapNum1` v

infixl 7 <.>, </>

-- | Gets direction angle (in radians) and produces a vector pointing in that direction
dirVec :: Double -> Vector
dirVec a = makeRel2 (cos a, sin a)

-- | Rotates the vector by the given angle
-- TODO: use rotateZaxis & multMatrix
rotateVec :: Double -> Vector -> Vector
rotateVec a v = makeRel2 (t `dotProduct` v, tp `dotProduct` v)
  where
    t = makeRel2 (cos a, -sin a)
    tp = makeRel2 (sin a, cos a)

-- | Converts an angle to the equivalent angle in (-pi,pi] range
fixAngle :: Double -> Double
fixAngle a  | a <= -pi  = fixAngle (a+2*pi)
            | a > pi    = fixAngle (a-2*pi)
            | otherwise = a

reflectVectorAgainst :: Vector -> Vector -> Vector
reflectVectorAgainst a b = c <.> 2 - a
  where
    au = unitVector a
    bu = unitVector b
    c = bu <.> (a `dotProduct` bu)

