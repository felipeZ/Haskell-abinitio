

module Science.QuantumChemistry.Indexes


import Data.Array.Repa  as R


-- R.fromIndex :: sh -> Int -> sh
-- R.toIndex :: sh -> sh -> Int 


toIndex4centers :: sh -> sh -> Int
toIndex4centers shape ix = R.fromIndex shape sh
 where sh = undefined

fromIndex4centers :: sh -> Int -> sh
fromIndex4centers shape i = R.fromIndex shape k
 where k =

sh4DtoList :: sh -> [Int]
sh4DtoList (Z :. i :. j :. k :. l) = [i,j,k,l]
