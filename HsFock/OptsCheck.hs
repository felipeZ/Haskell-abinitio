
-- Monadic command line checkers..
-- @2012,2013 Angel Alvarez Adhesive tape
-- @2012,2013 Felipe Zapata core SCF machinery 

module OptsCheck where

-- Main module needed imports
import Control.Monad(foldM,liftM,ap)
import Control.Monad.IO.Class
import Control.Monad.Trans.Either
import Data.List (find)
import Data.Maybe ( fromMaybe )
import Data.Either

import System.Environment
import System.Console.GetOpt
import System.Directory ( doesDirectoryExist, doesFileExist )
import System.FilePath


import Options


-- An EitherT container to store parsed opts from commandline or error messages
type OptsResult = EitherT String IO Options

-- An Opts filter runnng in the EitherT IO Stack
type OptsFilter = ( Options -> OptsResult )

-- A Policy describing a command line options with a checking filter
type OptsPolicy = OptDescr OptsFilter


-- =============================================== Options checking ==============================================
-- progOpts: getOpt args processor that also checks semantically getopt results or bang in front of the user
-- Upon checking with getopt this function gets a list of lambdas representing semantical checks
-- on every cmdline switch, as an example; we check for input file presence and the check that either data
-- director is a valid one and input file exists and is readable. Those checks are performed by
-- filename_check, check_datadir and check_input_file respectively. These funs are stored
-- in the Options structure that getopt uses.

--  We use a EitherT transformer to combine Either chaining with arbitrary IO actions needed during checking
-- ===============================================================================================================

progOpts :: [String] -> Options -> [OptsPolicy] -> OptsResult
progOpts args defaultOptions acceptedOptions =
   case getOpt RequireOrder acceptedOptions args of
--       (funs,[],[]) -> do
--           left "input file(s) missing"
      (funs,filenames,[]) -> do
          resultOfFuns <- foldl (>>=) (return defaultOptions) funs               -- Perform monadic checkings upon getOpt supplied functions
          foldM check_input_file resultOfFuns $ reverse filenames                          -- Now check if all the input files exist and are accesible
      (_,_,errs) -> do
          left ( concat errs )

--check supplied input files exist or bang if any is not.
check_input_file :: Options -> String -> OptsResult
check_input_file optsR@Options { optInput = files , optDataDir = dataDir } filename = do
    test <- liftIO $ filename_check dataDir filename
    case test of
         True -> return $ (optsR { optInput = filename : files })
         False -> left $ "input file "++ filename ++ " not readable"
    where
        filename_check :: Maybe FilePath -> FilePath -> IO Bool --check file with or without data directory
        filename_check (Just datadir) filename = doesFileExist $ combine datadir filename
        filename_check Nothing filename        = doesFileExist filename


