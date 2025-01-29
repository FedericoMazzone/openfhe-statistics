#!/bin/bash

# Check if the correct number of arguments are passed
if [ "$#" -ne 1 ]; then
    echo "Error: Please provide one argument: [ranking, ranking-tie, minimum, median, sorting]."
    exit 1
fi

# Assign the argument to a variable
algorithm=$1

# Define the corresponding executable based on the argument
tieParam=-1
case $algorithm in
    "ranking")
        executable="./build/ranking"
        tieParam=0
        ;;
    "ranking-tie")
        executable="./build/ranking"
        tieParam=1
        ;;
    "minimum")
        executable="./build/minimum"
        ;;
    "median")
        executable="./build/median"
        ;;
    "sorting")
        executable="./build/sorting"
        ;;
    *)
        echo "Error: Invalid argument. Please use one of the following: [ranking, ranking-tie, minimum, median, sorting]."
        exit 1
        ;;
esac

# Initialize output file
output_file="benchmark.out"

# if file is empty or does not exist, write the header
if [ ! -s "$output_file" ]; then
    echo "algorithm,thread,length,runtime(s),memory(MB)" > "$output_file"
fi

# Create a directory for logs if it does not exist
if [ ! -d "logs" ]; then
    mkdir logs
fi

# Iterate through the parameters from 2 to 6
for threadParam in 0 1; do
    for lengthParam in 8 16 32 64 128 256 512 1024 2048 4096 8192 16384; do

        if [ "$threadParam" -eq 1 ]; then
            threadString="single"
        else
            threadString="multi"
        fi
        
        echo "Running $algorithm on vector of size $lengthParam in $threadString-thread..."

        # Initialize a log file for the output of the program
        log_file="logs/${algorithm}-${threadString}thread-${lengthParam}.log"

        # If tieParam is defined, pass it as an argument
        if [ "$tieParam" -eq -1 ]; then
            command="$executable $lengthParam $threadParam"
        else
            command="$executable $lengthParam $tieParam $threadParam"
        fi
        echo $command
        # Start the program in the background, redirecting both stdout and stderr to the log file
        $command > "$log_file" 2>&1 &
        pid=$!
    
        # Initialize variables for tracking
        max_res=0
    
        # Monitor the RES (Resident Memory) every 0.5s
        while kill -0 "$pid" 2>/dev/null; do
            # Get the current RES in KB
            current_res=$(ps -o rss= -p "$pid" 2>/dev/null)
            
            # Ensure current_res is not empty or null
            current_res=${current_res:-0}
        
            # Update max_res if current_res is larger
            if [ "$current_res" -gt "$max_res" ]; then
                max_res="$current_res"
            fi
        
            # Sleep for 0.5 seconds before next check
            sleep 0.5
        done
    
        # Once the process ends, calculate max_res in MB
        final_max_res=$(awk "BEGIN {printf \"%.2f\", $max_res/1024}")

        # In the log file, search for the runtime of the program
        runtime=$(grep -oP 'Runtime: \K[0-9]+\.[0-9]+' "$log_file")
        
        # Append the result to the output file
        echo "$algorithm,$threadString,$lengthParam,$runtime,$final_max_res" | tee -a "$output_file"

    done
done

echo "Benchmarking completed. Results saved in $output_file."
