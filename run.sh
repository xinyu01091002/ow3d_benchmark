#!/bin/bash

# Environment variable for controlling filtered tasks
RUN_FILTERED_TASKS=${RUN_FILTERED_TASKS:-false}

# Log file
log_file="task_log.txt"

# CPU & Memory Threshold
CPU_THRESHOLD=95
MEMORY_THRESHOLD=80

# Delay time between tasks (in seconds)
TASK_DELAY=5

# Monitoring interval (in seconds)
RESOURCE_CHECK_INTERVAL=10

# Check completion interval (5 minutes)
COMPLETION_CHECK_INTERVAL=300  # 5*60 seconds

# Email Configuration
EMAIL="thooo2616@outlook.com"
EMAIL_SUBJECT="All Tasks Completed"
EMAIL_BODY="All tasks have been completed. CPU usage is below 1%."

# Check CPU and Memory usage
check_resources_usage() {
    while true; do
        cpu_usage=$(top -bn2 | grep "Cpu(s)" | sed -n '2p' | awk '{print 100 - $8}' | cut -d. -f1)
        memory_used_percent=$(free | awk '/Mem:/ {printf "%.0f", $3/$2 * 100}')

        echo "Current CPU usage: ${cpu_usage}% | Memory usage: ${memory_used_percent}%" | tee -a "$log_file"

        if [[ "$cpu_usage" -ge "$CPU_THRESHOLD" || "$memory_used_percent" -ge "$MEMORY_THRESHOLD" ]]; then
            echo "** Resource usage high (CPU: ${cpu_usage}%, Memory: ${memory_used_percent}%), waiting... **" | tee -a "$log_file"
            sleep $RESOURCE_CHECK_INTERVAL
        else
            echo "Resource usage is within limits (CPU: ${cpu_usage}%, Memory: ${memory_used_percent}%), proceeding with task creation." | tee -a "$log_file"
            break
        fi
    done
}

# Launch tasks
for dir in */ ; do
    dir_name=${dir%/}
    
    if [[ "$RUN_FILTERED_TASKS" == "true" ]]; then
        should_run=false
        for num in "${numbers_to_match[@]}"; do
            if [[ "$dir_name" =~ Akp_0$num ]]; then
                should_run=true
                break
            fi
        done
        if ! $should_run; then
            continue
        fi
    fi

    check_resources_usage

    echo "Waiting for ${TASK_DELAY} seconds before starting task in $dir_name..." | tee -a "$log_file"
    sleep $TASK_DELAY

    echo "Launching task in $dir_name" | tee -a "$log_file"
    task_start_time=$(date +"%Y-%m-%d %H:%M:%S")

    # Start task in screen
    screen -dmS "${dir_name}" bash -c "cd '$dir_name' && OceanWave3D_Tim; exec bash" &
    echo "Task started for $dir_name at $task_start_time" | tee -a "$log_file"
done

# Wait for all background tasks to finish launching
wait

# Monitor CPU usage periodically
echo "All tasks launched. Monitoring CPU usage every 5 minutes..." | tee -a "$log_file"

while true; do
    cpu_usage=$(top -bn2 | grep "Cpu(s)" | sed -n '2p' | awk '{print 100 - $8}' | cut -d. -f1)
    echo "$(date): CPU Usage: ${cpu_usage}%" | tee -a "$log_file"

    if [[ "$cpu_usage" -lt 1 ]]; then
        echo "CPU usage below 1%. Tasks are likely completed." | tee -a "$log_file"
        
        # Send email notification
        echo "$EMAIL_BODY" | mail -s "$EMAIL_SUBJECT" "$EMAIL"
        
        echo "Email notification sent to $EMAIL." | tee -a "$log_file"
        break
    fi

    echo "Tasks are still running. Checking again in 5 minutes..." | tee -a "$log_file"
    sleep $COMPLETION_CHECK_INTERVAL
done

# Final log entry
echo "All tasks completed. Total time elapsed: $SECONDS seconds." | tee -a "$log_file"
