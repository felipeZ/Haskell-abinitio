
-- The HaskellFock SCF Project 
-- @2013 Angel Alvarez, Felipe Zapata
-- Boys function calculation and related routines aimed at acurate and efficien 
-- calculation of Gaussian Type Orbitals (GTO)


-- Boys function calculation to multielectron Integrals
-- For a review about this approach, please refer to: 
--                 Multi-electron Integrals
--                 Simen Reine,1 Trygve Helgaker1 and Roland Lindh2
--                 WIREs Comput Mol Sci 2012, 2: 290–303 doi: 10.1002/wcms.78

-- [1] Kummer hypergeometric function using gaussian quadrature (adapted to use of external package)
-- For a complete description, please refer to: 
--                 Gauss quadrature approximations to hypergeometric and confluent hypergeometric functions
--                 Walter Gautschi
--                 Journal of Computational and Applied Mathematics 139 (2002) 173–187

-- [2] Intermediate range boys evaluation proposed on:
--                On the evaluation of Boys functions using downward recursion relations
--                B.A. Mamedov
--                Journal of Mathematical Chemistry Vol. 36, No. 3, July 2004 (© 2004)


-- Boys Function Fm(x) comparison against data from reference [2]
--         m = 8.0         x = 16.0        Eq(3) = 4.0230859250266e-7      Eq(6) = 4.0230859250266e-7      ours = 4.0230961756435804e-7
--         m = 15.0        x = 27.0        Eq(3) = 1.08359515555596e-11    Eq(6) = 1.08359515555596e-11    ours = 1.0836002458122441e-11
--         m = 20.0        x = 30.0        Eq(3) = 1.37585444267909e-3     Eq(6) = 1.37585444267909e-3     ours = 1.3758628962386526e-13
--         m = 25.0        x = 13.0        Eq(3) = 8.45734447905704e-8     Eq(6) = 8.45734447905704e-8     ours = 8.45737576279772e-8
--         m = 31.0        x = 34.0        Eq(3) = 2.90561943091301e-16    Eq(6) = 2.90561943091301e-16    ours = 2.9056447142297187e-16
--         m = 11.0        x = 38.0        Eq(3) = 4.04561442253925e-12    Eq(6) = 4.04561442253925e-12    ours = 4.0456286019975055e-12
--         m = 42.0        x = 32.0        Eq(3) = 5.02183610419087e-16    Eq(6) = 5.02183610419086e-16    ours = 5.021881853038279e-16
--         m = 75.0        x = 30.0        Eq(3) = 1.01429517438537e-15    Eq(6) = 1.01429517438537e-15    ours = 1.0143042539734162e-15
--         m = 100.0       x = 33.0        Eq(3) = 3.42689684943483e-17    Eq(6) = 3.42689684943483e-17    ours = 3.426930819254854e-17
--         m = 20.0        x = 1.4e-3      Eq(3) = 2.43577075309547e-2     Eq(6) = 2.43577075309547e-2     ours = 2.435770754086431e-2
--         m = 45.0        x = 6.4e-5      Eq(3) = 1.09883228385254e-2     Eq(6) = 1.09883228385254e-2     ours = 1.098832283873514e-2
--         m = 100.0       x = 2.6e-7      Eq(3) = 4.97512309732144e-3     Eq(6) = 4.97512309732144e-3     ours = 4.975123097321828e-3

module Science.QuantumChemistry.NumericalTools.Boys where

import Data.Foldable (mapM_)
import Data.Number.Erf
import Data.Vector.Unboxed as U
import Math.Gamma
import Prelude hiding (mapM_)

-- Internal modules
import Science.QuantumChemistry.GlobalTypes (VecUnbox)
import Science.QuantumChemistry.NumericalTools.PointsWeights 

-- Boys function
boysF :: Double -> Double -> Double
boysF 0 x
    | x > 0.0e-5 = 0.5 * sqrt (pi/x) * erf (sqrt x)
    | otherwise  = 1.0

boysF n x = k `seq` k / d
    where
        k = kummer a b z 0.1 -- Kummer 1F1
        a = n + 0.5
        b = n + 1.5
        z = -x
        d = (2 * n) + 1

-- Boys function Upward Recurrence as stated on [2]
-- boysF2 0 x = 
-- boysF2 m x = a * b
--     where
--         a = 1 / ( 2*x )
--         b = (2*m - 1) * boysF2 (m-1) x - (e ** (-x))
--         e = exp 1
-- Bous function backward recurrence as stated on [2]
-- boysF3 = bF
--     where
--         bF m x = a * b
--             where
--                 a = 1 / (2*m -1)
--                 b = 2 * x * bF (m+1) x - (e ** (-x))
--                 e = exp 1



-- Kummer's' "1F1" a.k.a M(a,b,z) Confluent Hypergeometric function as in [1]
-- Aproximation by Gaussian Quadrature method from 128 up to 1024 points of resolution
kummer :: Double -> Double -> Double -> Double -> Double
kummer a b z err = gammaFactor * integralPart
    where
        gammaFactor       = (gamma b) / (gamma a * gamma (b-a))
        integralPart      = integrator err fun 0 1
        fun               = (\t -> (e ** (z * t)) * (1-t) ** (b-a-1) * t ** (a-1))
        e                 = exp 1
        integrator err
                | err > 0.1   = nIntegrate128
                | err > 0.01  = nIntegrate256
                | err > 0.001 = nIntegrate512
                | otherwise   = nIntegrate1024
{- INLINE kummer-}

-- --  Test values from [2] for various proposed aproximations of the boys Function
-- --  Equations 3 and 6
-- testdata :: [(Double,Double,Double,Double)]
-- testdata = 
--     [
--          (  8,     16, 4.02308592502660e-07, 4.02308592502660e-07 )
--         ,( 15,     27, 1.08359515555596e-11, 1.08359515555596e-11 )
--         ,( 20,     30, 1.37585444267909e-03, 1.37585444267909e-03 )
--         ,( 25,     13, 8.45734447905704e-08, 8.45734447905704e-08 )
--         ,( 31,     34, 2.90561943091301e-16, 2.90561943091301e-16 )
--         ,( 11,     38, 4.04561442253925e-12, 4.04561442253925e-12 )
--         ,( 42,     32, 5.02183610419087e-16, 5.02183610419086e-16 )
--         ,( 75,     30, 1.01429517438537e-15, 1.01429517438537e-15 )
--         ,(100,     33, 3.42689684943483e-17, 3.42689684943483e-17 )
--         ,( 20, 1.4e-3, 2.43577075309547e-02, 2.43577075309547e-02 )
--         ,( 45, 6.4e-5, 1.09883228385254e-02, 1.09883228385254e-02 )
--         ,(100, 2.6e-7, 4.97512309732144e-03, 4.97512309732144e-03 )
--     ]

-- test :: IO ()
-- test = do
--     putStrLn "Boys Function Fm(x) comparison against data from reference [2]"
--     mapM_ report testdata
--     where
--         report :: (Double,Double,Double,Double) -> IO ()
--         report (a,b,c,d) = do
--             let ours = boysF a b
--             putStrLn $ "\tm = " ++ (show a) 
--                     ++ "\tx = " ++ (show b) 
--                     ++ "\tEq(3) = " ++ (show c) 
--                     ++ "\tEq(6) = " ++ (show d) 
--                     ++ "\tours = " ++ (show ours)
              
                  

baseCase :: (Double -> Double) -> VecUnbox -> VecUnbox -> Double
baseCase func points weights = U.sum $ U.zipWith (\x w -> w*(func x + func(-x))) points weights

nIntegrate128 :: (Double -> Double) -> Double -> Double -> Double
nIntegrate128 func a b = 0.5*(b-a) * (baseCase (\x -> func $ 0.5*((b-a)*x+b+a)) points128 weights128)

nIntegrate256 :: (Double -> Double) -> Double -> Double -> Double
nIntegrate256 func a b = 0.5*(b-a) * (baseCase (\x -> func $ 0.5*((b-a)*x+b+a)) points256 weights256)

nIntegrate512 :: (Double -> Double) -> Double -> Double -> Double
nIntegrate512 func a b = 0.5*(b-a) * (baseCase (\x -> func $ 0.5*((b-a)*x+b+a)) points512 weights512)

nIntegrate1024 :: (Double -> Double) -> Double -> Double -> Double
nIntegrate1024 func a b = 0.5*(b-a) * (baseCase (\x -> func $ 0.5*((b-a)*x+b+a)) points1024 weights1024)
