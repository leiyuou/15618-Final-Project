import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class addWeight {
    static TreeMap<Integer, Map<Integer, Integer>> nodes = new TreeMap<>();
    static int nodeCount = 0;
    static int edgeCount = 0;
    public static void main(String[] args) throws IOException {
        readFile("/Users/leiyuou/Desktop/15618/fakefile/src/roadNet-PA.txt");
        generateRandomWeight();
        writeToFile("/Users/leiyuou/Desktop/15618/fakefile/src/roadNet-PA-weighted_50000.txt");
    }

    public static void readFile(String path) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(path));
        String str;
        while ((str = in.readLine()) != null) {
            String[] tmp = str.split("\t");
            int node1 = Integer.valueOf(tmp[0]);
            int node2 = Integer.valueOf(tmp[1]);
            if (node1 > 49999 || node2 > 49999) {
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
                    int randomWeight = (int)(Math.random() * 1000) + 1;
                    subEntry.setValue(randomWeight);
                    nodes.get(subEntry.getKey()).put(entry.getKey(), randomWeight);
                }
                edgeCount++;
            }
        }
        nodeCount = nodes.size();

    }

    public static void writeToFile(String path) throws IOException {
        BufferedWriter out = new BufferedWriter(new FileWriter(path));
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
