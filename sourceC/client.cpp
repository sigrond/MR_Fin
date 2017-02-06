//#include <sys/types.h>
//#include <sys/socket.h>
//#include <sys/un.h>
#include<WinSock2.h>
#include <cstdio>
//#include <unistd.h>
#include <cstdlib>
int main(int argc, char ** argv)
{
	int sock;
	struct sockaddr server;
	sock = socket(AF_UNIX, SOCK_STREAM, 0);
	if (sock < 0) {
		perror("opening stream socket");
		exit(1);
	}
	//server.sun_family = AF_UNIX;
	//strcpy(server.sun_path, argv[1]);


	if (connect(sock, (struct sockaddr *) &server, sizeof(struct sockaddr)) < 0) {
		closesocket(sock);
		perror("Error while connecting to a stream socket, did you start a server?");
		exit(1);
	}
	int dummy =1;
	if ( send(sock, &dummy, 1*sizeof(int),0) < 0 )
		perror("Error writing to a stream socket");
	if (recv(sock, &dummy, 1*sizeof(int),0) < 0 )
		perror("Error reading stream stocket");
	closesocket(sock);
	return 0;
}
