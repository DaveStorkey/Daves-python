#!/bin/bash
#
#

while getopts f:g:v:M:m:Xe:o:C:K opt
do
  case $opt in
      f) moo_list_file=$OPTARG  ;;
      g) moo_list_file2=$OPTARG  ;;
      v) varnames=$OPTARG ;;
      M) maskfile=$OPTARG ;;
      m) masknames=$OPTARG ;;
      X) invert_mask="true" ;;
      e) end_file=$OPTARG  ;;
      o) file_out=$OPTARG  ;;
      C) chunksize=$OPTARG  ;;
      K) keep_files="true" ;;
  esac
done
shift `expr $OPTIND - 1`

if [[ -z "$moo_list_file" ]]; then
    echo "Error: must specify at least one moo_list file."
    exit 11
fi

if [[ -z "$varnames" ]]; then
    echo "Error: must specify at least one variable to process."
    exit 11
else
    varnames=$(echo $varnames | sed 's/,/ /g')
fi

if [[ -z "$file_out" ]]; then
    echo "Error: must specify the output file stem."
    exit 11
fi

mask_options=""
if [[ -n "$maskfile" ]]; then
    mask_options="-M $maskfile"
    if [[ -n "$masknames" ]]; then
	mask_options="$mask_options -m $(echo $masknames | sed 's/,/ /g')"
    fi
    if [[ "$invert_mask" == "true" ]]; then
        mask_options="$mask_options -X"	    
    fi
fi

if [[ -z "$chunksize" ]]; then chunksize=10; fi

if [[ -n "$end_file" ]]
then
    if [[ "$(echo $end_file | cut -c1-3)" == "moo" || "$(echo $end_file | cut -c1)" == ":" ]]
    then
       # (might) need to restore the file from MASS
       end_file_basename=$(basename $end_file)
       if [[ ! -e "$end_file_basename" ]]
       then
          echo "Restoring $end_file"
          moo get $end_file .
       fi
       end_file=$end_file_basename
    fi
fi
   
moofilelist=""
filelist=""
# the sed command removes any empty lines at the end of the file
moo_file_last=$(cat $moo_list_file | sed '/^[[:space:]]*$/d' | tail -n1)
moofilelist2=""
infilelist2=()
filelist2=""
# set to a non-empty value so the if-test further down doesn't get tripped if file2=""
moo_file_last2="elephant"  
if [[ -n "$moo_list_file2" ]];then
    for file in $(cat $moo_list_file2 | sed '/^[[:space:]]*$/d' )
    do
	infilelist2+=($file)
    done
    moo_file_last2=$(cat $moo_list_file2 | sed '/^[[:space:]]*$/d' | tail -n1)
fi

let countchunk=0
let counttotal=0
append_option=""
for file in $(cat $moo_list_file | sed '/^[[:space:]]*$/d' )
do
    if [[ ${#infilelist2[@]} -ge $((counttotal+1)) ]];then
	file2=${infilelist2[$counttotal]}
    else
	file2=""
    fi
    ((counttotal++))
    ((countchunk++))
    if [[ ! -e $(basename $file) ]]
    then 
        moofilelist="$moofilelist $file"
    fi
    filelist="$filelist $(basename $file)"
    if [[ -n "$file2" ]];then
        if [[ ! -e $(basename $file2) ]]
        then 
            moofilelist2="$moofilelist2 $file2"
        fi
        filelist2="$filelist2 $(basename $file2)"
    fi
	
    if [[ $countchunk == $chunksize || "$file" == $moo_file_last || "$file2" == "$moo_file_last2" ]]
    then
        if [[ -n "$moofilelist" ]]
	then
           moo get $moofilelist .
        fi
        if [[ -n "$moofilelist2" ]]
	then
           moo get $moofilelist2 .
        fi
        filelist2_option=""
        if [[ -n "$filelist2" ]];then filelist2_option="-j $filelist2" ; fi
        end_file_option=""
        if [[ -n "$end_file" ]];then  end_file_option="-e $end_file" ; fi
        echo "filelist : $filelist"
        echo "filelist2 : $filelist2"
        calc_rms_series.py -i $filelist $filelist2_option -v $varnames $mask_options -o $file_out $end_file_option $append_option
        # overwrite output files on first call and append on subsequent calls to calc_rms_series.py
        append_option="-A"
        let countchunk=1
        moofilelist=""
        if [[ "$keep_files" != "true" ]];then
            filelist_to_delete=$(echo $filelist | rev | cut -d" " -f2- | rev)
            if [[ -n "$filelist2" ]];then
	        filelist_to_delete="$filelist_to_delete $(echo $filelist2 | rev | cut -d" " -f2- | rev)"
	    fi
	    rm $filelist_to_delete
	fi
        filelist=$(echo $filelist | rev | cut -d" " -f1 | rev)
        if [[ -n "$filelist2" ]];then
            moofilelist2=""
	    filelist2=$(echo $filelist2 | rev | cut -d" " -f1 | rev)
	    if [[ "$filelist2" == "$moo_file_last2" ]];then filelist2="";fi
	fi
    fi
done

exit 0

