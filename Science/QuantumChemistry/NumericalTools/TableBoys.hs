
-- The HaskellFock SCF Project 
-- @2013 Angel Alvarez, Felipe Zapata
-- Boys function calculation and related routines aimed at acurate and efficien 
-- calculation of Gaussian Type Orbitals (GTO)


-- Boys function calculation to multielectron Integrals
-- For a review about this approach, please refer to: 
--                 Multi-electron Integrals
--                 Simen Reine,1 Trygve Helgaker1 and Roland Lindh2
--                 WIREs Comput Mol Sci 2012, 2: 290–303 doi: 10.1002/wcms.78

module Science.QuantumChemistry.NumericalTools.TableBoys (
                                                         Boys
                                                        ,boysTaylor
                                                        ,generateGrid
                                                         ) where

import Control.Applicative ((<$>),(<*>))
import Data.Foldable (mapM_)
import Data.Map.Lazy as M
import Prelude hiding (mapM_)
import Text.Printf 

-- Internal modules
import Science.QuantumChemistry.GlobalTypes (VecUnbox)
import Science.QuantumChemistry.NumericalTools.PointsWeights
import Science.QuantumChemistry.NumericalTools.Boys (boysF) 

-- Data Types
data Boys = Boys Double Double deriving (Eq,Ord) -- ^represents the function f(m,x)


{- |
The boys function is cached on the fly using a lazy map that store the points of a grid.
The point f(m,x) is calculated as using a Taylor series expansion around xk,
f(m,x) = f(m,xk + dx) = f(m,xk) -f(m+1,xk)*dx + 0.5*f(m+2,xk)(dx)^2 -1/6 * f(m+3,xk)(dx)^2 +
         1/24 * f(m+4,xk)(dx)^2
  
-}
boysTaylor :: Map Boys Double -> Double -> Double -> Double
boysTaylor grid m x = (fun m xk) - (fun (m+1) xk)*dx + 0.5*(fun (m+2) xk)*dx^2 - (1/6)*(fun (m+3) xk)*dx^3 + (1/24)*(fun (m+4) xk)*dx^4 - (1/120)*(fun (m+5) xk)*dx^5
  where delta = 0.1
        dx    = x-xk
        xk    = (calcClosestPoint x )* delta
        fun   = calcBoys grid  
        calcClosestPoint r = fromIntegral  $ floor $ x /delta
       

generateGrid :: Int -> Double -> Map Boys Double
generateGrid mMax delta = M.fromList $ zip keys val 
  where keys = Boys <$> [fromIntegral x | x <- [0..mMax]] <*> [delta*fromIntegral(i) | i <- [0..100000]]
        val  = fmap (\(Boys m x) -> boysF m x) keys
         
calcBoys :: Map Boys Double -> Double -> Double -> Double
calcBoys grid m x =
   let msg = printf "The requested Boys function is outside the chosen grid: boysF %f %f" m x
   in case M.lookup (Boys m x) grid of
           Just v  -> v
           Nothing -> error msg

