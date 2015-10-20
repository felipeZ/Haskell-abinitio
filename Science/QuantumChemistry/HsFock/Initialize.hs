{-# LANGUAGE RecordWildCards #-}

module Science.QuantumChemistry.HsFock.Initialize where 

-- =============================> Standard and third party libraries <=================
import Control.Exception (throw)
import qualified Data.Map as M
import Text.Printf

-- =================> Internal Modules <======================
import Science.QuantumChemistry.Error (HSFOCKException(..))
import Science.QuantumChemistry.GlobalTypes (AtomData(..), CGF(..), HSFOCK(..), ZNumber, atomLabel2Charge)
import Science.QuantumChemistry.ParsecTools.ParseXYZ (parseFileXYZ)


-- ============================================================

-- | using the provided basis set and the molecular coordinates initialize the atom Data type
initializeAtoms :: HSFOCK -> (String -> IO () ) -> IO [AtomData]
initializeAtoms HSFOCK{..} logger =
 do  atomsXYZ        <- parseFileXYZ xyz logger 
     let basisSetMap =  createBasisMap atomsXYZ basis 
     return $ createAtomData atomsXYZ basisSetMap 


createAtomData :: [(String,[Double])] -> M.Map String [CGF] -> [AtomData]
createAtomData atomsXYZ basisSetMap = map create atomsXYZ
 where create (l,xyz) = AtomData xyz (lookupAtom l basisSetMap)  (lookupAtom l atomLabel2Charge)
       lookupAtom k m   = case M.lookup k m of
                          Nothing -> throw  KeyError       
                          Just v  ->  v 

createBasisMap :: [(String,[Double])] -> String -> M.Map String [CGF] 
createBasisMap atomsXYZ basis = undefined


