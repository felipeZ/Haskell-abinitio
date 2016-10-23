{-# LANGUAGE LambdaCase #-}

{-|
Module: Science.QuantumChemistry.BasisSet.FetchBasis
Description: Configure the basis set for the simulation
Copyright: @2016 Felipe Zapata
-}

module Science.QuantumChemistry.BasisSet.FetchBasis where

-- ==================> Standard and third party libraries <=================
import Control.Monad (unless, when)
import qualified Data.ByteString.Char8 as B
import Data.Char (toLower)
import qualified Data.Map.Strict  as M
import Data.Maybe (fromMaybe)
import Data.Serialize (decode)
import System.Directory (doesDirectoryExist, getDirectoryContents)
import System.Environment (getEnv)
import System.Exit (exitFailure)
import System.FilePath.Posix ((</>))
import Text.Printf

-- =================> Internal Modules <======================
import Science.QuantumChemistry.GlobalTypes (CGF(..), Funtype(..))
import Science.QuantumChemistry.ParsecTools.ParserBasis 
import Science.QuantumChemistry.BasisSet.NormalizeBasis (normGlobal)
-- ==========================<>=============================

{-| For a given basis, creates a map with Atom symbols as
keywords and a contracted Gaussian
functions as values.-}
createBasisMap :: [String] -> String -> IO (M.Map String [CGF])
createBasisMap xs basis =   do
  allElems  <- fetchBasis basis
  let mapa = createElemMap allElems
  return $ foldr (go mapa) M.empty xs

   where go elemMap x acc =
            case M.lookup x elemMap of
                 Just (Atom _ gs) -> M.insert x (createCGF gs) acc
                 Nothing          -> acc 


createElemMap ::[ Element] -> M.Map String Element
createElemMap = foldr (\e@(Atom label _) acc -> M.insert (B.unpack label) e acc) M.empty

{-| Using The Coefficients and Exponents create a basis Set, represented as a list
of contracted Gauss functions -}
createCGF :: [GaussShape] -> [CGF]
createCGF = concatMap fun 
  where  cgf xs s = [normGlobal (CGF xs s)]
         fun      = \case 
               SP  xs -> let cs1 = fmap (\(c1, _, e) -> (c1,e)) xs
                             cs2 = fmap (\(_, c2, e) -> (c2,e)) xs  
                         in fmap normGlobal [CGF cs1 S, CGF cs2 Px, CGF cs2 Py, CGF cs2 Pz]
               S0   xs -> [normGlobal (CGF xs S)]
               P    xs -> fmap normGlobal [CGF xs Px, CGF xs Py, CGF xs Pz]
               D    xs -> fmap normGlobal [CGF xs Dxx, CGF xs Dxy, CGF xs Dxz, CGF xs Dyy, CGF xs Dyz, CGF xs Dzz]


-- | Reads and Decodes a basis set store as Bytestring
fetchBasis :: String -> IO [Element]
fetchBasis basis = do
     pathBasis <- defaultPathBasis
     names     <- getDirectoryContents pathBasis
     let basisFiles = filter (`notElem` [".", ".."]) names
         realName   = translateBasisname basis
     when (realName `notElem` basisFiles) $  printOnExit $ printf "Unkown basis set:%s\n" basis
     xs <- B.readFile (pathBasis </> realName)
     case decode xs of
          Left msg -> printOnExit msg 
          Right rs -> return rs


-- | Default path where the basis sets are stored
defaultPathBasis :: IO FilePath
defaultPathBasis = do
  hsFockPath <- getEnv "HSFOCK"
  boolDir    <- doesDirectoryExist hsFockPath
  unless boolDir $ printOnExit "HsFock installation directory was not found, please add the following environmental variable to your profile\nHSFOCK=<path/to/the/HSFOCK/"
  return (hsFockPath </> "data/basis")

     
{-| Because Symbols like '+*' are not allowed as part of names, the basis set containing
these symbols are stored in files where `+` is replace by `m` and `*` by `s`.
For example the basis 6-31+G** is stored in the file `6-31mGss.basis`-}
translateBasisname :: String -> String
translateBasisname = (++ ".basis") . fmap (replace . toLower)
 where rs        = [('+','m'),('*','s')]
       replace x = fromMaybe x (lookup x rs) 

printOnExit :: String -> IO a
printOnExit s = putStrLn s *> exitFailure 
