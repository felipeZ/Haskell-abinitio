
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
                                                        ,generateGridBoys
                                                         ) where

import Control.Applicative ((<$>),(<*>))
import Data.Foldable (mapM_)
import Data.Map.Lazy as M
import Prelude hiding (mapM_)
import Text.Printf 

-- Internal modules
import Science.QuantumChemistry.GlobalTypes (VecUnbox)
import Science.QuantumChemistry.NumericalTools.PointsWeights
import Science.QuantumChemistry.NumericalTools.Boys (asymptBoysF,boysF) 

-- Data Types
data Boys = Boys Double Double deriving (Eq,Ord) -- ^represents the function f(m,x)


{- |
We can expand the Boys function as a 6 terms Taylor series around xk,
f(m,x) = f(m,xk + dx) = f(m,xk) -f(m+1,xk)*dx + 0.5*f(m+2,xk)(dx)^2 -1/6 * f(m+3,xk)(dx)^2 +
         1/24 * f(m+4,xk)(dx)^2
Using the Taylor series we can generate a grid of xk points [0.1, 0.2, ..xk,.. xmax] in mmax + 6 dimesions, where xmax and mmax are the maximum value of x and m, respectively. Basically, the idea is to populate the nodes of the grid calling the expensive function f and the value of the function between nodes are calculated using the above Taylor series expansion.
-}

boysTaylor :: Map Boys Double -> Double -> Double -> Double
boysTaylor grid m x = (fun m xk) - (fun (m+1) xk)*dx + 0.5*(fun (m+2) xk)*dx^2 - (1/6)*(fun (m+3) xk)*dx^3 + (1/24)*(fun (m+4) xk)*dx^4 - (1/120)*(fun (m+5) xk)*dx^5
  where delta = 0.1
        dx    = x-xk
        xk    = (calcClosestPoint x )* delta
        fun   = calcBoys grid  
        calcClosestPoint r = fromIntegral  $ floor $ x /delta
       

{- |
Generate a grid where each node is a thunk of the expensive function boysF. Now every time that we need to calculate boysF we call funTaylor instead of f. Notice that's this is a virtual grid in the sense that every time that we call funTaylor m x, only the points involved in the Taylor series are evaluated while the rest of the grid remains unevaluated.

-}
generateGridBoys ::  Double -> Map Boys Double
generateGridBoys delta = M.fromList $ zip keys val 
  where keys = Boys <$> [fromIntegral x | x <- [0..mMax]] <*> [delta*fromIntegral i | i <- [0..xMax]]
        val  = fmap (\(Boys m x) -> boysF m x) keys
        mMax = 14
        xMax = 1000
         
calcBoys :: Map Boys Double -> Double -> Double -> Double
calcBoys grid m x
  | x > 100  = asymptBoysF m x
  | otherwise =               
        let msg = printf "The requested Boys function is outside the chosen grid: boysF %f %f" m x
        in case M.lookup (Boys m x) grid of
                Just v  -> v
                Nothing -> error msg

