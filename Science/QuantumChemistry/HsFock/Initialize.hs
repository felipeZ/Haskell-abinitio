{-# LANGUAGE RecordWildCards #-}

module Science.QuantumChemistry.HsFock.Initialize where 

-- =============================> Standard and third party libraries <===============================


-- =================> Internal Modules <======================
import Science.QuantumChemistry.GlobalTypes
import Science.QuantumChemistry.ParsecTools.ParseXYZ (parseFileXYZ)

-- ============================================================

initializeAtoms :: HSFOCK -> IO [AtomData]
initializeAtoms HSFOCK{..} = undefined
  -- do
  -- atomsXYZ <- parseFileXYZ xyz
  
