

module Science.QuantumChemistry.NumericalTools.Utils (
             fac
            ,facOdd
            ,normGlobal
             )where


import Control.Arrow ((&&&))
import Control.Applicative ((<$>),(<*>))
import Data.Foldable (sum)
import Prelude hiding (sum)

-- Internal modules 
import Science.QuantumChemistry.GlobalTypes
import Science.QuantumChemistry.NumericalTools.LinearAlgebra as LA
import Science.QuantumChemistry.Integrals.IntegralsEvaluation (fac, facOdd, sijContracted)

-- | see http://www.ccl.net/cca/documents/basis-sets/basis.html 
normGlobal :: CGF -> CGF
normGlobal cgf@(CGF ps l) =  CGF newPrimitives l
  where (cs,es)       = unzip . getPrimitives $ cgf
        xs            = normCoeff cgf
        cgfN          =  CGF (zip es xs) l
        r0            = [0,0,0]
        normaG        = sqrt . recip $  sijContracted r0 r0 cgfN cgfN 
        newCs         = map (*normaG) xs                
        newPrimitives = zip newCs es 


{- |The norm of each contracted is given by the following equation
    N = sqrt $ ((2l -1)!! (2m-1)!! (2n-1)!!)/(4*expo)^(l+m+n)  * (pi/(2*e))**1.5
    where expo is the exponential factor of the contracted -}
normCoeff :: CGF -> [Double]
normCoeff b1 = zipWith fun cs es 
  where (cs,es) = unzip . getPrimitives $ b1
        fun = \c e -> (c/) . sqrt $ (ang e) * (pi/(2*e))**1.5
        ang x = prod / (4*x)^(sum indexes)
        prod = product $ fmap (\k -> facOdd (2*k -1)) indexes
        shell = getfunTyp  b1
        indexes = fmap (LA.map2val mapLAngular) $ zip (repeat shell) [0..2]

-- normCoeff :: CGF -> CGF
-- normCoeff b1 = b1 { getPrimitives = newPrimitives}
--   where xs = getPrimitives b1
--         newPrimitives = fmap ((uncurry fun ) &&& (snd) ) xs
--         fun = \c e -> (c/) . sqrt $ (ang e) * (pi/(2*e))**1.5
--         ang x = prod / (4*x)^(sum indexes)
--         prod = product $ fmap (\k -> facOdd (2*k -1)) indexes
--         shell = getfunTyp  b1
--         indexes = fmap (LA.map2val mapLAngular) $ zip (repeat shell) [0..2]


-- normCoeff :: Funtype -> [GaussPrimitive] -> [Double]
-- normCoeff orbType gs = map fun gs
--   where fun = \(c,e) -> (c/) . sqrt $ (ang e) * (pi/(2*e))**1.5
--         ang x = prod / (4*x)^(sum indexes)
--         prod = product $ fmap (\k -> facOdd (2*k -1)) indexes
--         indexes = fmap (LA.map2val mapLAngular) $ zip (repeat orbType) [0..2]



