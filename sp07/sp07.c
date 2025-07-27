#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#define LARGE_DATA_SIZE 5000000
#define SLEEP_TIME 2

void print_timestamp(int rank, const char* message) {
    static double start_time = 0.0;
    if (start_time == 0.0) {
        start_time = MPI_Wtime();
    }
    double current_time = MPI_Wtime() - start_time;
    printf("[%.2fs] Rank %d: %s\n", current_time, rank, message);
    fflush(stdout);
}

void heavy_computation(int rank, const char* phase) {
    print_timestamp(rank, "💻 Starting heavy computation...");
    
    // 重い計算をシミュレート（実際の計算処理）
    long result = 0;
    for (int i = 0; i < 100000000; i++) {
        result += i % 1000;
    }
    
    char msg[256];
    sprintf(msg, "✅ Computation completed (result: %ld) - %s", result % 10000, phase);
    print_timestamp(rank, msg);
}

void test_blocking_communication(int rank, int size) {
    printf("\n");
    printf("========================================\n");
    printf("🔴 BLOCKING COMMUNICATION TEST (MPI_Send)\n");
    printf("========================================\n");
    fflush(stdout);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) {
        int *data = (int*)malloc(LARGE_DATA_SIZE * sizeof(int));
        for (int i = 0; i < LARGE_DATA_SIZE; i++) {
            data[i] = i;
        }
        
        print_timestamp(0, "📤 About to call MPI_Send...");
        double start = MPI_Wtime();
        
        // ブロッキング送信
        MPI_Send(data, LARGE_DATA_SIZE, MPI_INT, 1, 100, MPI_COMM_WORLD);
        
        double send_time = MPI_Wtime() - start;
        char time_msg[256];
        sprintf(time_msg, "📤 MPI_Send COMPLETED after %.3f seconds", send_time);
        print_timestamp(0, time_msg);
        
        print_timestamp(0, "⚠️  Notice: I was BLOCKED until send completed!");
        
        // 送信完了後に計算
        heavy_computation(0, "AFTER blocking send");
        
        free(data);
        
    } else if (rank == 1) {
        print_timestamp(1, "😴 Simulating slow receiver (sleeping 3 seconds)...");
        sleep(3);  // 受信を意図的に遅らせる
        
        int *recv_data = (int*)malloc(LARGE_DATA_SIZE * sizeof(int));
        
        print_timestamp(1, "📥 About to receive data...");
        MPI_Recv(recv_data, LARGE_DATA_SIZE, MPI_INT, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        print_timestamp(1, "📥 Data received successfully!");
        
        free(recv_data);
    }
}

void test_nonblocking_communication(int rank, int size) {
    printf("\n");
    printf("========================================\n");
    printf("🟢 NON-BLOCKING COMMUNICATION TEST (MPI_Isend)\n");
    printf("========================================\n");
    fflush(stdout);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) {
        int *data = (int*)malloc(LARGE_DATA_SIZE * sizeof(int));
        for (int i = 0; i < LARGE_DATA_SIZE; i++) {
            data[i] = i + 1000;
        }
        
        print_timestamp(0, "📤 About to call MPI_Isend...");
        double start = MPI_Wtime();
        
        MPI_Request request;
        // ノンブロッキング送信
        MPI_Isend(data, LARGE_DATA_SIZE, MPI_INT, 1, 200, MPI_COMM_WORLD, &request);
        
        double isend_time = MPI_Wtime() - start;
        char time_msg[256];
        sprintf(time_msg, "📤 MPI_Isend RETURNED immediately (%.6f seconds)", isend_time);
        print_timestamp(0, time_msg);
        
        print_timestamp(0, "🚀 I can do other work while sending!");
        
        // 送信中に他の作業ができる
        heavy_computation(0, "DURING non-blocking send");
        
        print_timestamp(0, "⏳ Now checking if send is complete...");
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        double total_time = MPI_Wtime() - start;
        
        sprintf(time_msg, "✅ Send actually completed after %.3f seconds total", total_time);
        print_timestamp(0, time_msg);
        
        free(data);
        
    } else if (rank == 1) {
        print_timestamp(1, "😴 Simulating slow receiver (sleeping 3 seconds)...");
        sleep(3);  // 受信を意図的に遅らせる
        
        int *recv_data = (int*)malloc(LARGE_DATA_SIZE * sizeof(int));
        
        print_timestamp(1, "📥 About to receive data...");
        MPI_Recv(recv_data, LARGE_DATA_SIZE, MPI_INT, 0, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        print_timestamp(1, "📥 Data received successfully!");
        
        free(recv_data);
    }
}

void show_summary(int rank) {
    if (rank == 0) {
        printf("\n");
        printf("========================================\n");
        printf("📊 SUMMARY: KEY DIFFERENCES\n");
        printf("========================================\n");
        printf("\n");
        printf("🔴 MPI_Send (BLOCKING):\n");
        printf("   ├─ Function does NOT return until send is safe\n");
        printf("   ├─ CPU waits and cannot do other work\n");
        printf("   ├─ Simpler to use (no need to check completion)\n");
        printf("   └─ May waste CPU cycles waiting\n");
        printf("\n");
        printf("🟢 MPI_Isend (NON-BLOCKING):\n");
        printf("   ├─ Function returns IMMEDIATELY\n");
        printf("   ├─ CPU can do computation while communication happens\n");
        printf("   ├─ Must use MPI_Wait/MPI_Test to check completion\n");
        printf("   └─ Better performance through computation/communication overlap\n");
        printf("\n");
        printf("💡 USE NON-BLOCKING WHEN:\n");
        printf("   ├─ You have computation to do while waiting\n");
        printf("   ├─ Sending large amounts of data\n");
        printf("   ├─ Need maximum performance\n");
        printf("   └─ Want to overlap multiple communications\n");
        printf("\n");
        printf("========================================\n");
        fflush(stdout);
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (size < 2) {
        if (rank == 0) {
            printf("❌ Error: This program requires at least 2 processes\n");
            printf("Usage: mpirun -np 2 ./program\n");
        }
        MPI_Finalize();
        return 1;
    }
    
    if (rank == 0) {
        printf("🚀 MPI BLOCKING vs NON-BLOCKING COMMUNICATION DEMO\n");
        printf("Data size: %d integers (%.2f MB)\n", LARGE_DATA_SIZE, 
               (LARGE_DATA_SIZE * sizeof(int)) / (1024.0 * 1024.0));
        printf("Processes: %d\n", size);
        printf("\nWatch the timestamps to see the difference!\n");
    }
    
    test_blocking_communication(rank, size);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("\n⏱️  Waiting 2 seconds between tests...\n");
        sleep(2);
    }
    test_nonblocking_communication(rank, size);
    
    MPI_Barrier(MPI_COMM_WORLD);
    show_summary(rank);
    
    MPI_Finalize();
    return 0;
}
