module Primitives (
  Arc(..),
  isArc, isSegment,
  arcFromPt, arcToPt,
  arcNormals,
  arcAngles
) where

import Data.SG (plusDir, perpendicular2, unitVector, fromPt, dotProduct)
import Geometry (Point, Vector, dirVec, (<.>), fixAngle)

data Arc  = A { center :: Point, radius :: Double, fromA :: Double, toA :: Double }
          | L { linePt1 :: Point, linePt2 :: Point }
          deriving Show

isArc, isSegment :: Arc -> Bool
isArc (A _ _ _ _) = True
isArc _ = False
isSegment (L _ _) = True
isSegment _ = False

arcFromPt, arcToPt :: Arc -> Point
arcFromPt arc | isArc arc     = (center arc) `plusDir` ((dirVec $ fromA arc) <.> (radius arc))
              | isSegment arc = linePt1 arc
arcToPt   arc | isArc arc     = (center arc) `plusDir` ((dirVec $ toA arc) <.> (radius arc))
              | isSegment arc = linePt2 arc

arcNormals :: Arc -> ((Point, Vector), (Point, Vector))
arcNormals arc | isArc arc = ((p1,v1), (p2, v2))
               | otherwise = ((p1,lv), (p2, lv))
  where
    p1 = arcFromPt arc
    p2 = arcToPt arc
    v1 = (perpendicular2 $ dirVec $ fromA arc) <.> signum a
    v2 = (perpendicular2 $ dirVec $ toA arc) <.> signum a
    a = fixAngle (toA arc - fromA arc)
    lv = unitVector $ p2 `fromPt` p1

-- debugging function
arcAngles :: Arc -> Arc -> Double
arcAngles a1 a2 = 180 * (acos $ (unitVector a2n) `dotProduct` (unitVector a1n)) / pi
  where
    a1n = snd $ snd $ arcNormals a1
    a2n = snd $ fst $ arcNormals a2

