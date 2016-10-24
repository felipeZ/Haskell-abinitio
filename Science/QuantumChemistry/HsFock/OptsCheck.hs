{-# LANGUAGE  RecordWildCards #-}

{-|
Module: Science.QuantumChemistry.HsFock.OptsCheck
Description: Command line checkers
Copyright: @2016 Felipe Zapata
-}

module Science.QuantumChemistry.HsFock.OptsCheck (checkOpts) where

-- =============================> Standard and third party libraries <===============================
import Control.Monad (unless, when)
import Data.Char (toLower)
import System.Directory ( doesDirectoryExist, doesFileExist, listDirectory )
import System.Exit (exitFailure)
import Text.Printf

-- =================> Internal Modules <=====================
import Science.QuantumChemistry.GlobalTypes (HSFOCK(..))

-- ==============================<>==============================

-- | Sanity check of hsfock arguments
checkOpts :: HSFOCK -> IO ()
checkOpts opts@HSFOCK{..}= do
  unless  scf       $ printExit "We only Serve HartreeFock SCF (for the moment)"
  when (null basis) $ printExit "You must provide a basis set, like: sto-3g, etc."
  existBasis        <- doesBasisExists basis
  unless existBasis $ printExit "Unknown basis set"
  when (charge < 0) $ printExit "Charge must be a natural number"
  when (multi <= 0) $ printExit "Multiplicity must be 1 (singlet), 2(doblet), etc."
  --todo check that basis set exist
  [xyzBool,outBool] <- mapM doesFileExist [xyz, outFile]
  unless  xyzBool . printExit $ printf "Coordinates File %s does not exist" xyz 
  when (outBool)  . printExit $ printf "file %s is going to be overwritten\n" outFile 

   where printExit s = putStrLn s *> exitFailure 

-- | Sanity check of basis set configuration
checkOpts opts@BasisConfig{..}= do
   fileBool <- doesFileExist basisPath
   when fileBool ( printf "The file:%s\n already exists, it is going to overwritten\n" basisPath)

-- | Check if the user provided basis exits
doesBasisExists :: String -> IO Bool
doesBasisExists bs = listBasis >>= \xs -> return (basis `elem` xs) 
  where basis = fmap toLower bs 

-- | List available basis set
listBasis :: IO [String]
listBasis = do
  xs <- listDirectory "data/basis"
  return (fmap (takeWhile (/= '.')) xs)
