{-# LANGUAGE RecordWildCards #-}

{-|
Module: Science.QuantumChemistry.HsFock.Initialize
Description: Initialize basis set for the atoms
Copyright: @2013 Felipe Zapata, Angel Alvarez
           @2016 Felipe Zapata
-}

module Science.QuantumChemistry.HsFock.Initialize where 

-- =============================> Standard and third party libraries <=================
import Control.Arrow (second)
import Control.Exception (throw)
import Data.List (nub)
import Data.Maybe (fromMaybe)
import qualified Data.Map.Strict as M
import Text.Printf

-- =================> Internal Modules <======================
import Science.QuantumChemistry.BasisSet.FetchBasis (createBasisMap)
import Science.QuantumChemistry.Error (HSFOCKException(..))
import Science.QuantumChemistry.GlobalTypes (AtomData(..), CGF(..), HSFOCK(..), ZNumber, atomLabel2Charge)
import Science.QuantumChemistry.ParsecTools.ParseXYZ (parseFileXYZ)
import Science.QuantumChemistry.Units (angstrom2AU)

-- ============================================================

-- | Using the provided basis set and the molecular coordinates initialize the atom Data type
initializeAtoms :: HSFOCK -> (String -> IO () ) -> IO [AtomData]
initializeAtoms HSFOCK{..} logger =
 do  atomsXYZ        <- (fmap (second (fmap angstrom2AU))) <$> parseFileXYZ xyz logger
     let uniqueElems = nub $ map fst atomsXYZ
     basisSetMap <- createBasisMap uniqueElems basis
     return $ createAtomData atomsXYZ basisSetMap 

-- | Create Atom types
createAtomData :: [(String,[Double])] -> M.Map String [CGF] -> [AtomData]
createAtomData atomsXYZ basisSetMap = map create atomsXYZ
 where create (l,xyz) = AtomData xyz (lookupAtom l basisSetMap)  (lookupAtom l atomLabel2Charge)
       lookupAtom k m = fromMaybe (throw KeyError) (M.lookup k m)




