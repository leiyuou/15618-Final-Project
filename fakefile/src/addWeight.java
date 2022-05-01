import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class addWeight {
    static final int NODE = 70000;
    static final int MAX = 200;
    static TreeMap<Integer, Map<Integer, Integer>> nodes = new TreeMap<>();
    static int nodeCount = 0;
    static int edgeCount = 0;
    public static void main(String[] args) throws IOException {
        readFile("C:/Users/xinqi/Desktop/15618/15618-Final-Project/fakefile/src/roadNet-PA.txt");
        generateRandomWeight();
        writeToFile(String.format("C:/Users/xinqi/Desktop/15618/15618-Final-Project/inputs/roadNet-PA-%d.txt", NODE));
    }

    public static void readFile(String path) throws IOException {
        FileReader fr = new FileReader(path);
        BufferedReader in = new BufferedReader(fr);
        String str;
        while ((str = in.readLine()) != null) {
            String[] tmp = str.split("\t");
            int node1 = Integer.valueOf(tmp[0]);
            int node2 = Integer.valueOf(tmp[1]);
            if (node1 > NODE-1 || node2 > NODE-1 || (Math.abs(node1-node2)<2)) {
                continue;
            }
            if (!nodes.containsKey(node1)) {
                nodes.put(node1, new HashMap<>());
            }
            Map<Integer, Integer> adjNodes = nodes.get(node1);
            adjNodes.put(node2, 0);
        }
        in.close();
    }

    public static void generateRandomWeight() {
        for (Map.Entry<Integer, Map<Integer, Integer>> entry : nodes.entrySet()) {
            Map<Integer, Integer> map = entry.getValue();
            for (Map.Entry<Integer, Integer> subEntry : map.entrySet()) {
                if (subEntry.getValue() == 0) {
                    int randomWeight = (int)(Math.random() * MAX) + 1;
                    subEntry.setValue(randomWeight);
                    nodes.get(subEntry.getKey()).put(entry.getKey(), randomWeight);
                }
                edgeCount++;
            }
            int u = entry.getKey();
            if(u-1>-1) {
                map.put(u - 1, (NODE / 2) * MAX);
                edgeCount++;
            }
            if(u+1<NODE) {
                map.put(u + 1, (NODE / 2) * MAX);
                edgeCount++;
            }
        }
        nodeCount = NODE;
    }

    public static void writeToFile(String path) throws IOException {
        File f = new File(path);
        f.createNewFile();
        BufferedWriter out = new BufferedWriter(new FileWriter(f));
        out.write(String.format("%d\t%d\n", nodeCount, edgeCount));
        for (Map.Entry<Integer, Map<Integer, Integer>> entry : nodes.entrySet()) {
            Map<Integer, Integer> map = entry.getValue();
            for (Map.Entry<Integer, Integer> subEntry : map.entrySet()) {
                out.write(String.format("%d\t%d\t%d\n", entry.getKey(), subEntry.getKey(), subEntry.getValue()));
            }
        }
        out.flush();
        out.close();
    }
}
