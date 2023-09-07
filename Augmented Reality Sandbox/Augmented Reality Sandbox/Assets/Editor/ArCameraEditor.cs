using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(ArCamera))]
public class ArCameraEditor : Editor
{
    string[] cameraDeviceNames;

    ArCamera arCamera;

    void Awake()
    {
        arCamera = GameObject.FindObjectOfType<ArCamera>();

        cameraDeviceNames = new string[WebCamTexture.devices.Length];

        for (int i = 0; i < WebCamTexture.devices.Length; i++)
        {
            cameraDeviceNames[i] = WebCamTexture.devices[i].name;
        }
    }

    void OnSceneGUI()
    {
        //Handles.ScaleSlider(1, Vector3.zero, Vector3.up, Quaternion.identity, 1, 0f);
    }
 
    public override void OnInspectorGUI()
    {
        DrawDefaultInspector();

        EditorGUILayout.Space();

        GUI.changed = false;

        arCamera.webCamIndex = EditorGUILayout.Popup("Camera Device", arCamera.webCamIndex, cameraDeviceNames, EditorStyles.popup);

        arCamera.usePrerecordedVideo = EditorGUILayout.Toggle("Use Prerecorder Video", arCamera.usePrerecordedVideo);
        arCamera.folderPath = EditorGUILayout.TextField("Source Video", arCamera.folderPath);


        if (GUILayout.Button("Select Folder"))
            arCamera.folderPath = EditorUtility.SaveFolderPanel("Save textures to directory", "C:\\Users\\Aytek Aman\\Dropbox\\Research\\Benchmark", "");

        if (GUI.changed)
            EditorUtility.SetDirty(arCamera);

        if (GUILayout.Button("Perturb"))
        {
            arCamera.transform.position = arCamera.testPosition + Random.insideUnitSphere * 0.05f;
            arCamera.transform.rotation = Quaternion.Euler(arCamera.testRotation + new Vector3(Random.Range(0f,1f), Random.Range(0f,1f), Random.Range(0f,1f)) * 5);
        }


        //EditorGUILayout.Space();
        //GUI.changed = false;
        //if (GUILayout.Button("Process Markers"))
        //{
        //    //arCamera.Init();
        //    //arCamera.CreateChildCameras();
        //    arCamera.CombineMarkers();
        //    arCamera.ProcessMarkers();
        //    arCamera.InitCells();
        //    arCamera.CalculateVisibility();
        //    //arCamera.CalculateReliability();
        //}

        //if (GUI.changed)
        //{
        //    EditorUtility.SetDirty(arCamera);
        //}
        
        //GUILayout.Button("Cache visibility");
        //GUILayout.Button("Process illumination");

        //GUI.changed = false;
        //arCamera.cellSize =  EditorGUILayout.FloatField("Cell Size", arCamera.cellSize);

        //if (GUI.changed)
        //    EditorUtility.SetDirty(arCamera);
    }
}