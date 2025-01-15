import pandas as pd
import numpy as np
import argparse
import re
import json
import os

##
# Return a dictionary
# @dtype Type of dictionary to return [acq, mod, label]
# @return dictionary in hashable form 
def getDictInfo(dtype,dicts):
    if dtype == "acq":
        info = os.path.join(dicts,"acq_dict.json")
    elif dtype == "mod":
        info = os.path.join(dicts,"mod_dict.json")
    elif dtype == "label":
        info = os.path.join(dicts,"label_dict.json")
    
    with open(info) as json_file:
        info_dict = json.load(json_file)

    return info_dict

##
# Return 0 instead of nan for a given value
# @val Any variable of any type
# @return Actual value if variable is not nan, if nan returns 0
def getVal(val):
  if pd.isnull(val) == False:
    return float(val)
  else:
    return 0


##
# Deal with uneven magnetic field strength numbers
# @param strength Numeric representation of the reported magnetic field strength
# @return correctedStrength
def correctMagneticFieldStrength(strength):
    # Make sure the strength is numeric
    strength = float(strength)

    if strength == 1.5 or strength == 3.0:
        return strength

    # Get the difference between reported strength and both options
    delta15 = np.abs(1.5 - strength)
    delta30 = np.abs(3.0 - strength)
    # How to generalize this for other T? - added to TODO list
    if delta15 < delta30:
        return 1.5
    elif delta30 < delta15:
        return 3.0
    else:
        print("Unable to correct field strength:", strength)
        return np.nan


##
# Remove any 2D images using the whole summary df
# @param summaryDf A Pandas dataframe containing a CuBIDS summary.tsv file data
# @return summaryDf Same as the input but modified to have MergeInto column values of 0 for 2D scans
def remove2DScans(summaryDf):
    summaryDf.loc[
        summaryDf["Dim4Size"] <= 1.0 and summaryDf["Dim3Size"] <= 1.0, "MergeInto"
    ] = 0
    return summaryDf


def filterDWI(df, paramGroup, idx):
   if paramGroup["Modality"] == "dwi":
      fieldStrength = str(correctMagneticFieldStrength(paramGroup["MagneticFieldStrength"])) + "T"
      # Drop anything with less than 7 volumes
      if paramGroup["NumVolumes"] < 7:
        df.loc[idx, "MergeInto"] = 0
      # Add number of volumes to filename
      newName = createNewName(paramGroup, paramGroup["ProtocolName"], fieldStrength, "dwi")
      df.loc[idx, "RenameEntitySet"] = newName
        

##
# Identify DTI map types: ADC, TRACEW, and MDDW
#   ADC:
#   TRACEW:
#   MDDW:
def identifyDTISubtypes(paramGroup):
    pass

##
# Compare a given scan's parameters to a pre-determined group of parameters for specific sequences
# @input is a dataframe containing pre-determined group of sequences with specific parameters
# @paramGroup is a single row in summary dataframe
# @flex_tr is the flex range for repetition time
# @flex_te is the flex range for echo time
# @return output Name(s) of protocols a given parameter group matches 
def comparison(input, paramGroup, flex_tr, flex_te):
  
  mfs = getVal(paramGroup["MagneticFieldStrength"])
  mfs = float(round(mfs,1))   
  fieldStrength = str(mfs)+"T"
  orient = paramGroup["ImageOrientation"]
  te = getVal(paramGroup["EchoTime"])
  tr = getVal(paramGroup["RepetitionTime"])
  flAng = getVal(paramGroup["FlipAngle"])
  perSam = getVal(paramGroup["PercentSampling"])
  mb_acc = getVal(paramGroup["MultibandAccelerationFactor"])
  pxb = getVal(paramGroup["PixelBandwidth"])
  output = input[(fieldStrength == input.MagneticFieldStrength)&
                             # Filter for fixed parameters
                            ((orient == input.ImageOrientation)|
                            (pd.isnull(orient) == True))&
                            ((flAng == input.FlipAngle)|
                            (flAng == input.FlipAngle_alt)|
                            (flAng == 0))&
                            ((perSam == input.PercentSampling)|
                            (perSam == input.PercentSampling_alt)|
                            (perSam == 0))&
                            (((pxb == input.PixelBandwidth)|
                              (pxb == input.PixelBandwidth_alt)|
                              (pxb == 0)))&
                            (((mb_acc == input.MultibandAccelerationFactor)|
                              (mb_acc == input.MultibandAccelerationFactor_alt)|
                              (mb_acc == 0)))&
                            # Filter for minimal tolerance parameters
                            (((tr <= (input.RepetitionTime + flex_tr))&
                            (tr >= (input.RepetitionTime - flex_tr)))|
                            (tr == input.RepetitionTime_alt))&
                            (((te <= (input.EchoTime + flex_te))&
                            (te >= (input.EchoTime - flex_te)))|
                            (te == input.EchoTime_alt))]["ProtocolName"].unique()
  return output

##
# Get the general acquisition type for a given parameter group (scan)
# @paramGroup A single row from the summary.tsv file
# @return prefix Acquisition type from the acqusition dictionary
def getGeneralCategory(paramGroup):
  key = paramGroup["EntitySet"]
  prefix = None
  if "acquisition" in key:
    prefix = key.split("acquisition-")[1].split("_")[0]
  return prefix

##
# Get suffix
# @paramGroup A single row from the summary.tsv file
# @grouping predefined scan groups
# @matches Group that a scan falls into after doing comparison
# @fieldStrength magnetic field strength
# @return suffix for a scan [T1w, T2w, FLAIR, unknown]
def getSuffix(paramGroup, grouping, matches, fieldStrength,dicts=False):
  keyGroup = paramGroup["EntitySet"]
  suffixes = []

  if not dicts:
    for s in matches:
      val = grouping[(grouping.ProtocolName == s)&
                          (fieldStrength == grouping.MagneticFieldStrength)]["group"].item()
      suffixes.append(val)     
      if len(set(suffixes))==1:
        suffix = "_suffix-"+suffixes[0]
      else:
        suffix = "_suffix-"+keyGroup.split("suffix-")[1]
  else: 
    if "unknown" in keyGroup:
      suffix = "_suffix-unknown"
      if matches:          
        # Use the label dictionary to get a suffix
        label_dict = getDictInfo("label", dicts)
        match = matches.lower()
        suffix = "_suffix-unknown" if not [k for k in label_dict if(k in match or k in match)] else "_suffix-%s" % label_dict[[k for k in label_dict if(k in match or k in match)][0]]
        
      # try coarse match 
      if "unknown" in suffix:
        suffix = "_suffix-%s" % CoarseMatch(paramGroup, grouping, fieldStrength)
    else:
      suffix = "_suffix-" + keyGroup.split("suffix-")[1]

  if "w1" in keyGroup or "w2" in keyGroup:
    ending = str(re.findall(r'w1|w2', keyGroup)[0].split("w")[1])
    suffix = suffix + ending
  return suffix


##
# Coarsely define an unknown scan as T1w, T2w, or FLAIR based only on TR, TE and Magnetic Field Strength
# @paramGroup A single row from the summary.tsv file
# @grouping predefined scan groups
def CoarseMatch(paramGroup, grouping, fieldStrength):
  tr = round(paramGroup["RepetitionTime"],1)
  te = round(paramGroup["EchoTime"],1)
  coarse_match = grouping[(grouping["MagneticFieldStrength"]==fieldStrength)&
                          (grouping["RepetitionTime"]>=(tr-0.05))&
                          (grouping["RepetitionTime"]<=(tr+0.05))&
                          (grouping["EchoTime"]>=(te-0.005))&
                          (grouping["EchoTime"]<=(te+0.05))]["group"].unique()
  if len(set(coarse_match)) == 1:
    return coarse_match[0]
  else:
    return "unknown"


def createNewName(paramGroup, pName, fieldStrength, suffix):
  mod = paramGroup["Modality"]
  keyGroup = paramGroup["EntitySet"]


  if mod == "anat":
    newGroupName = [
        None if not "sample" in keyGroup else "sample-fetal_",
        "acquisition-" if not pName else "acquisition-%s" % pName,
        fieldStrength,
        None if not paramGroup["DeviceSerialNumber"] else "ScannerId%s" % str(paramGroup["DeviceSerialNumber"]),
        None if not paramGroup["Dim3Size"] else "Slices%d" % int(paramGroup["Dim3Size"]),
        None if not paramGroup["VoxelSizeDim1"] else "Voxel%smm" % str(round(float(paramGroup["VoxelSizeDim1"]),3)),
        None if not "ceagent" in keyGroup else "_%s" % "ceagent-gad",
        "_datatype-anat",
        None if not "reconstruction" in keyGroup else "_%s" % re.findall("reconstruction-[A-Z0-9]+",keyGroup)[0], 
        None if not "echo" in keyGroup else "_%s" % re.findall(r'echo-[0-9]+',keyGroup)[0],
        # None if not "part" in keyGroup else "_%s" % re.findall(r'part-[mag|phase|real|imaginary]',keyGroup)[0],
        None if not "mt" in keyGroup else "_mt-on",
        "_run-"+re.findall(r'run-[0-9]+',paramGroup["EntitySet"])[0].split("-")[1].lstrip('0').zfill(2),
        suffix
      ]
  
  elif mod == "dwi":
    newGroupName = [
        None if not "sample" in keyGroup else "sample-fetal_",
        "acquisition-%s" % fieldStrength,
        None if not paramGroup["DeviceSerialNumber"] else "ScannerId%s" % str(paramGroup["DeviceSerialNumber"]),
        None if not paramGroup["NumVolumes"] else "Volumes%d" % int(paramGroup["NumVolumes"]),
        None if not "ceagent" in keyGroup else "_%s" % "ceagent-gad",
        "_datatype-dwi",
        None if not "reconstruction" in keyGroup else "_%s" % re.findall("reconstruction-[A-Z0-9]+",keyGroup)[0], 
        None if not "mt" in keyGroup else "_mt-on",
        "_run-"+re.findall(r'run-[0-9]+',paramGroup["EntitySet"])[0].split("-")[1].lstrip('0').zfill(2),
        "_suffix-dwi"
      ]
    
  protocol = "".join(filter(bool, newGroupName)).replace("3.0T","3T").replace(".","p")
  return protocol
##
# Identify scans using predefined protocol groups and parameters to determine a standard or variant type.
# @paramGroup A single row from the summary.tsv file
# @grouping Table with predefined acquisition types
# @return newGroupName A string representing a new name for files in the specified group OR np.nan
#   File name: acquisition-[Protocol][Standard/Variant](details)_[bids_fields]_run-(runNumber)_[suffix]
def identify_standard_variant(paramGroup,grouping,dicts):
 
  keyGroup = paramGroup["EntitySet"]
  mfs = correctMagneticFieldStrength(paramGroup["MagneticFieldStrength"])
  dim3 = getVal(paramGroup["Dim3Size"])
  fieldStrength = str(mfs)+"T"
  
  standard_scans = []
  variant_scans = []
  protocol = None
  
  if dim3 > 10 and mfs:
   

    if "MPRAGEStandardized" not in keyGroup:
      standards = comparison(grouping,paramGroup,0.05,0.01)
      # print("Standards:",*standards)
      
      
      if len(set(standards)) == 1:
        standard = grouping[(grouping.ProtocolName == standards[0])&
                            (fieldStrength == grouping.MagneticFieldStrength)]["ShortName"].item()
        standard_scans.append(standard)
        protocol = "MPRAGEVariant" if "MPR" in standard else standard+"StandardizedVariant"
        # Get suffix
        suffix = "_suffix-"+grouping[(grouping.ProtocolName == standards[0])&
                            (fieldStrength == grouping.MagneticFieldStrength)]["group"].item()

      elif len(set(standards)) > 1:
        protocol = "MPRAGEVariant" if all(['MPR' in item for item in standards]) else getGeneralCategory(paramGroup)
        # Check if suffix the same for all groups
        suffix = getSuffix(paramGroup, grouping,standards,fieldStrength)
      
      
      else:
        variants = comparison(grouping,paramGroup,0.5,0.05)
        # print("Variants:",*variants)
        
        
        if len(set(variants)) == 1:
          variant = grouping[(grouping.ProtocolName == variants[0])&
                            (fieldStrength == grouping.MagneticFieldStrength)]["ShortName"].item()
          variant_scans.append(variant)
          protocol = variant+"Variant"
          suffix = "_suffix-"+grouping[(grouping.ProtocolName==variants[0])&(grouping.MagneticFieldStrength==fieldStrength)]["group"].item()

        
        elif len(set(variants)) > 1:
          protocol = "MPRAGEVariant" if all(['MPR' in item for item in variants]) else getGeneralCategory(paramGroup)
          # Check if suffix the same for all groups
          suffix = getSuffix(paramGroup, grouping,variants,fieldStrength)

        else:
          protocol = getGeneralCategory(paramGroup)
          suffix = getSuffix(paramGroup, grouping, protocol, fieldStrength, dicts)
      
      protocol = createNewName(paramGroup, protocol, fieldStrength, suffix)
      return protocol



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", 
        "--summary-file", 
        help="Full path to the run01_summary.tsv file."
    )
    parser.add_argument(
        "-fn", 
        "--output-file", 
        help="Output filename.",
        required=False
    )
    parser.add_argument(
        "-g",
        "--groups",
        help="Full path to the tsv file containing predefined scan groups.",
    )
    parser.add_argument(
        "-d",
        "--dictionary",
        help="Full path to the driectory containing the label dictionary.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print visual feedback",
    )
    parser.add_argument(
        "-m",
        "--modality",
        help="Only process specific modalities either anat, dwi or both. Defaults to only process anat data.",
        default="anat",
        type=str,
        nargs='+'
    )

    args = parser.parse_args()
    fn = args.summary_file
    gs = args.groups
    dicts = args.dictionary

    # Read in a summary.tsv file
    df = pd.read_csv(fn, sep="\t")
    preDefinedGroups = pd.read_csv(gs, sep="\t")
    print("input summary table", df.shape)

    # For my sanity:
    df = df[df["EchoTime"].notna()]
    print("input summary tale after filtering for EchoTime", df.shape)

    # For each parameter group in the summary dataframe
    for idx, group in df.iterrows():
      newName = None
      kg = group["EntitySet"]

      toIgnore = re.findall(r'w\d+|e\d+|Standardized|Variant|Scanner',kg)
      if group["Modality"] != "dwi":
        if ('anat' in kg or 'unknown' in kg) and not toIgnore:
          newName = identify_standard_variant(group, preDefinedGroups, dicts)
              
        # Make manual notes
        if newName:
          
          if "unknown" not in newName:
            if  "unknown" in kg:
              # Label groups where the coarse search was used and returned a suffix [T1w,T2w, or FLAIR] based on field strength, TR and TE
              if "acquisition-1" in newName or "acquisition-3" in newName:
                df.loc[idx, "Notes"] = "coarse mfs/tr/te match"
              else:
                # Label groups where the suffix was determined by the acquisition label (i.e. MTS usually T1w, SPC usually T2w)
                if "Variant" not in newName:
                  df.loc[idx, "Notes"] = "heuristic label match"
          else:
            # We don't want CuBIDS to act on the files that are still failing to be classified and remain "unknown"
            newName = None
          
        

        # Fix issue where files with "echo" in the name get overwritten:
        if "echo" in kg and newName:
          echoNum = newName.split("echo-")[1].split("_")[0]
          newEnd = "w" + str(echoNum)
          newName = newName.replace("w", newEnd)

        df.loc[idx, "RenameEntitySet"] = newName
      
      elif group["Modality"] == "dwi" and "dwi" in args.modality:
        filterDWI(df, group, idx)


    # Visual feedback: print non-nan param group names
    if args.verbose:
      print(df[df["RenameEntitySet"].notna()]["RenameEntitySet"].values)

    fnOut = fn.replace("_summary","_summary_edited") if not args.output_file else args.output_file
    df.to_csv(fnOut, sep="\t", index=False)


if __name__ == "__main__":
    main()