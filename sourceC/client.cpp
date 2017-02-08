//#include <sys/types.h>
//#include <sys/socket.h>
//#include <sys/un.h>
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <WinSock2.h>
#include <ws2tcpip.h>
#include <iphlpapi.h>
#include <cstdio>
//#include <unistd.h>
#include <cstdlib>
#define DEFAULT_PORT "27015"
#pragma comment (lib, "Ws2_32.lib")
#pragma comment (lib, "Mswsock.lib")
#pragma comment (lib, "AdvApi32.lib")
int main(int argc, char ** argv)
{
	int sock;
	struct sockaddr server;
	//sock = socket(AF_UNIX, SOCK_STREAM, 0);
	//if (sock < 0) {
	//	perror("opening stream socket");
	//	exit(1);
	//}
	//server.sun_family = AF_UNIX;
	//strcpy(server.sun_path, argv[1]);

	WSADATA wsaData;
	int iResult;
	// Initialize Winsock
	iResult = WSAStartup(MAKEWORD(2, 2), &wsaData);
	if (iResult != 0)
	{
		printf("WSAStartup failed: %d\n", iResult);
		return 1;
	}

	struct addrinfo *result = NULL, *ptr = NULL, hints;

	ZeroMemory(&hints, sizeof(hints));
	hints.ai_family = AF_INET;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_protocol = IPPROTO_TCP;
	hints.ai_flags = AI_PASSIVE;

	iResult = getaddrinfo(NULL, DEFAULT_PORT, &hints, &result);
	if (iResult != 0) {
		printf("getaddrinfo failed with error: %d\n", iResult);
		WSACleanup();
		return 1;
	}

	sock = socket(result->ai_family, result->ai_socktype, result->ai_protocol);
	if (sock == INVALID_SOCKET)
	{
		printf("Error at socket(): %ld\n", WSAGetLastError());
		freeaddrinfo(result);
		WSACleanup();
		return 1;
	}


	iResult = connect(sock, ptr->ai_addr, (int)ptr->ai_addrlen);
	if (iResult == SOCKET_ERROR) {
		closesocket(sock);
		sock = INVALID_SOCKET;
		return 1;
	}

	char dummy[sizeof(int)];
	if ( send(sock, dummy, 1*sizeof(int),0) < 0 )
		perror("Error writing to a stream socket");
	if (recv(sock, dummy, 1*sizeof(int),0) < 0 )
		perror("Error reading stream stocket");
	closesocket(sock);
	freeaddrinfo(result);
	return 0;
}
