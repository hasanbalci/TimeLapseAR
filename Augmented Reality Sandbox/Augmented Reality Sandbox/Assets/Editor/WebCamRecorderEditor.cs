using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(WebCamRecorder))]
public class WebCamRecorderEditor : Editor
{
    string[] cameraDeviceNames;

    WebCamRecorder webCamRecorder;

    void Awake()
    {
        webCamRecorder = GameObject.FindObjectOfType<WebCamRecorder>();

        cameraDeviceNames = new string[WebCamTexture.devices.Length];

        for (int i = 0; i < WebCamTexture.devices.Length; i++)
        {
            cameraDeviceNames[i] = WebCamTexture.devices[i].name;
        }
    }

    public override void OnInspectorGUI()
    {
        DrawDefaultInspector();

        GUI.changed = false;

        webCamRecorder.webCamIndex = EditorGUILayout.Popup("Camera Device", webCamRecorder.webCamIndex, cameraDeviceNames, EditorStyles.popup);

        if (GUILayout.Button("Select Folder"))
            webCamRecorder.folderPath = EditorUtility.SaveFolderPanel("Select Folder", "C:\\Users\\Aytek Aman\\Dropbox\\Research\\Benchmark", "Recording");

        if (GUI.changed)
            EditorUtility.SetDirty(webCamRecorder);
        
        //EditorGUILayout.Space();

        //GUI.changed = false;

        //arCamera.webCamIndex = EditorGUILayout.Popup("Camera Device", arCamera.webCamIndex, cameraDeviceNames, EditorStyles.popup);

        //if (GUI.changed)
        //    EditorUtility.SetDirty(arCamera);
    }
}