module Math (
  expm1,
  coshm1,
  squareEqSolutions
) where

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


-- | Solve a x^2 + b x + c == 0
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

