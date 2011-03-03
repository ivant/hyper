module Main where

import Geometry ((<.>), fixAngle, dirVec, Point, rotateVec)
import Primitives (arcNormals, center, radius, isArc, toA, fromA, linePt1, linePt2, Arc(..), arcFromPt, arcToPt)
import Hyperbolic (reflectThrough, arcFrom, reflectPtThrough)

import Data.SG (plusDir, Pair(..), iso, getX, getY, Point2'(..), unitVector)
import Graphics.Rendering.Cairo (Render, lineTo, moveTo, stroke, save, restore, setSourceRGB, setDash, arc, arcNegative, scale, translate, setLineWidth, paint, withImageSurface, Format(..), renderWith, surfaceWriteToPNG)
import Text.Printf (printf)
import Data.Maybe (fromJust)

type Ngon = [Arc]

setupPNG :: Int -> Render ()
setupPNG sz = do
    let sc = 1.3
    scale (fromIntegral sz / (2*sc)) (-fromIntegral sz / (2*sc))
    translate (sc) (-sc)
    setLineWidth $ ((2*sc) / fromIntegral sz)

reflectNgonThrough :: Ngon -> Arc -> Ngon
reflectNgonThrough p a = map (`reflectThrough` a) p

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

    drawNgon (0,0,1) ngon
    let ngons = do m <- ngon
                   return $ ngon `reflectNgonThrough` m

    mapM_ (drawNgon (0,0.5,0.7)) ngons
    
    return ()

  where
    drawNgon c p = mapM_ (drawArc c) p

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

    drawArc (r,g,b) a | isArc a = do
                          drawArcConstruction a

                          save
                          setDash [] 0
                          setSourceRGB r g b
                          let arcF = if (fixAngle $ toA a - fromA a) > 0 then arc else arcNegative
                          arcF (getX (center a)) (getY (center a)) (radius a) (fromA a) (toA a)
                          stroke
                          restore
                          --drawNormals a
                      | otherwise = drawStraight (r,g,b) (linePt1 a) (linePt2 a)

ngon = makeNgon 5 (pi/2) (\r -> Point2 (0, 0)) (pi/4+pi/2)

makeNgon :: Int -> Double -> (Double -> Point) -> Double -> [Arc]
makeNgon ngonSides ngonAngle startPF startAngle = take ngonSides $ iterate (fromJust . nextArc ngonAngle len) firstArc
  where
    firstArc = fromJust $ arcFrom startP startDir
    --firstArc = fromJust $ arcFrom (Point2 (0,0.1)) (dirVec 0 <.> len)
    startP = startPF r
    startDir = dirVec startAngle <.> len

    r = let t = tan (beta/2) in sqrt ((1-t)/(1+t))
    nextArc alpha len arc = arcFrom (fst n) dir
      where
        n = snd $ arcNormals arc
        dir = (rotateVec alpha (unitVector $ snd n)) <.> len

    beta = 2*pi/(fromIntegral ngonSides)
    len = acosh (1+2*(cos beta))



main = do
  let sz = 600
  withImageSurface FormatRGB24 sz sz $ \s -> do
    renderWith s $ (setupPNG sz >> drawHyper)
    surfaceWriteToPNG s "hyper.png"

  return ()
