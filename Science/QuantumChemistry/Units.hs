{-|
Module      : Science.QuantumChemistry.Units
Description : Back and for unit transformation functions
-}
module Science.QuantumChemistry.Units where 


-- =========================<>=================================
-- | Transform back and form between atomic unit of 
au2Angstrom, angstrom2AU :: Double -> Double

au2Angstrom = (*0.529177)
angstrom2AU = (/0.529177)
