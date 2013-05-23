import gazebo.WorldStatistics;
import gazebo.Time;

import java.io.BufferReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.IOException
import java.io.PrintStream;

class Server {

	boolean paused = false;
	int iterations = 222;


	void publishWorldStats() {
		
		WorldStatistics.Builder worldStats = WorldStatistics.newBuilder();
		Time.Builder simTime = Time.newBuilder();
		Time.Builder pauseTime = Time.newBuilder();
		Time.Builder realTime = Time.newBuilder();

		simTime.setSec(20);
		simTime.setNsec(21);

		pauseTime.setSec(30);
		pauseTime.setNsec(31);

		realTime.setSec(40);
		realTime.setNsec(41);

		worldStats.setSim_time(simTime);	
		worldStats.setPause_time(pauseTime);	
		worldStats.setReal_time(realTime);	
		worldStats.setPaused(paused);
		worldStats.setIteration(iteration);

	}

	void serve() {

		while(true)
		{

		}
	}
		
	void connect() {
		String address = "127.0.0.1";
		int port = 11345;
		Socket s;

		try {
			s = new Socket(address, port);
		} catch (IOException e) {
			e.printStackTrace();
		}	
		
		BufferedOutputStream bis = new BufferedOutputStream(s.getOutputStream());

		try {
			bis.write()
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}

public static void main(String[] args) throws Exception {

	Server s = new Server();
	s.serve();
}
