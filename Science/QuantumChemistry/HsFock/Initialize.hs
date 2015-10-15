{-# LANGUAGE RecordWildCards #-}

module Science.QuantumChemistry.HsFock.Initialize where 

-- =============================> Standard and third party libraries <=================
import Control.Exception (throwIO)
import qualified Data.Map as M


-- =================> Internal Modules <======================
import Science.QuantumChemistry.Error (HSFOCKException(..))
import Science.QuantumChemistry.GlobalTypes (AtomData, CGF(..), HSFOCK(..),ZNumber)
import Science.QuantumChemistry.ParsecTools.ParseXYZ (parseFileXYZ)

-- ============================================================

-- | using the provided basis set and the molecular coordinates initialize the atom Data type
initializeAtoms :: HSFOCK -> (String -> IO () ) -> IO [AtomData]
initializeAtoms HSFOCK{..} logger =
 do  atomsXYZ        <- parseFileXYZ xyz logger 
     let basisSetMap =  createBasisMap atomsXYZ basis 
     createAtomData atomsXYZ basisSetMap logger 


createAtomData :: [(String,[Double])] -> M.Map String [CGF] -> (String -> IO ()) -> IO [AtomData]
createAtomData atoms basisSetMap logger  = undefined
 -- mapM create atoms
 -- where create (l,xyz) = AtomData xyz
 --      lookupIO k m    =
 --               case M.lookup k m of
 --                    Nothing -> do let err = printf "unkown atom:%s\n" k
 --                                  logger (printf ) *> throwIO KeyError       
 --                    Just v  -> return v 

createBasisMap :: [(String,[Double])] -> String -> M.Map String [CGF] 
createBasisMap = undefined


