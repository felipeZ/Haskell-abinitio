
-- Monadic command line checkers..
-- HaskellFock SCF Project 
-- @2012 Angel Alvarez
-- @2012 Felipe Zapata

module OptsCheck where

-- Main module needed imports
-- import Control.Monad
import Control.Monad(foldM,liftM,ap)
import Data.Maybe ( fromMaybe )
import Data.Either
import System.Environment
import System.Console.GetOpt
import System.Directory ( doesDirectoryExist, doesFileExist )
import System.FilePath.Posix

-- Record for storing cmdline options
data Options = Options
 { optDump        :: Bool
 , optVerbose     :: Bool
 , optShowVersion :: Bool
 , optOutput      :: Maybe FilePath
 , optInput       :: Maybe FilePath
 , optDataDir     :: Maybe FilePath
 } deriving Show

-- Monadic Check functions for arguments, we live in the IO Monad for now mostly to allow us
-- manage stuff like file IO and exception handling...
fun_check_data_dir :: String 
    -> Either String Options 
    -> IO ( Either String Options )
fun_check_data_dir dir opts = do
    case opts of
         Left _msg -> return opts
         Right optsR -> do
             test <- doesDirectoryExist dir
             case test of
                  True -> return $ Right ( optsR { optDataDir =  Just dir } )
                  False -> return $ Left ( "Data directory " ++ dir ++ " does not exist" )

-- check user wants verbosity level increased
fun_check_verbosity :: Either String Options 
    -> IO( Either String Options )
fun_check_verbosity opts = do
    case opts of
         Left _msg -> return opts
         Right optsR -> return $ Right ( optsR { optVerbose = True } )

-- Who we are?, sort of alien outbreak?
fun_check_version :: Either String Options 
    -> IO ( Either String Options )
fun_check_version opts = do
    case opts of
         Left _msg -> return opts
         Right optsR ->  return $ Right ( optsR { optShowVersion = True } )

-- dump either options or errors as we get passthrought
fun_check_dump_options :: Either String Options 
    -> IO ( Either String Options)
fun_check_dump_options opts = do
    case opts of
         Left msg -> do
             putStrLn $ "Options got error: " ++ msg ++ "\n"
             return opts
         Right optsR -> do
             putStrLn $ "Options dumping selected record: " ++ show optsR ++ "\n"
             return $ Right ( optsR { optDump =True } )

-- help message
fun_check_help :: Either String Options 
    -> IO ( Either String Options)
fun_check_help opts = do
    return $ Left "Command line help requested"

--check file with or without data directory
filename_check :: Maybe FilePath 
    -> FilePath 
    -> IO Bool
filename_check (Just datadir) filename = do
    doesFileExist $ combine datadir filename
filename_check Nothing filename = do
    doesFileExist filename

--check basis file exists or use default "input.dat"
check_input_file :: [String] 
    -> Either String Options 
    -> IO (Either String Options)
check_input_file _ (Left msg) = do
    return $ Left msg
check_input_file [] _ = do
    return $ Left "input file missing"
check_input_file (filename:_) (Right optsR) = do
    test <- filename_check (optDataDir optsR) filename
    case test of
         True -> return $ Right (optsR { optInput = Just filename})
         False -> return $ Left "input file not readable"

-- =============================================== Options checking ==============================================
-- progOpts: getOpt args processor that also checks semantically getopt results or bang in front of the user
-- Upon checking with getopt this function gets a list of lambdas representing semantical checks
-- on every switch, as an example; we check for input file presence and the check that either data
-- director is a valid one and input file exists and is readable. Those checks are performed by
-- filename_check, fun_check_datadir and fun_check_input_file respectively. These funs are stored
-- in the Options structure that getopt uses.

progOpts :: [String] 
    -> Options 
    -> [OptDescr (Either String Options -> IO (Either String Options))] 
    -> String
    -> IO (Either String Options)
progOpts argv defaultOptions acceptedOptions header =
   case getOpt RequireOrder acceptedOptions argv of
      (funs,filenames,[]) -> do
          -- Inverse fold over the checking funs... (A nasty try to get a Monadic Either experience....)
          result <- foldM (flip id) ( Right defaultOptions) funs
          -- Now check if the input file exist and is accesible
          check_input_file filenames result
      (_,_,errs) -> do
          return $ Left ( concat errs ) 