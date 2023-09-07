using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System;

public class WebCamRecorder : MonoBehaviour
{
    public enum WebCamResolution
    {
        Low,
        Medium,
        High
    }

    //[HideInInspector]
    public int webCamIndex;

    int width, height;

    public WebCamResolution webCamResolution = WebCamResolution.High;

    bool isRecording = false;

    WebCamTexture webCamTexture;

    List<Color32[]> frames;

    public string folderPath;

    byte[] buffer;

    int frameCount = 0;

    int bufferSize;

    float recordingStartTime, recordingEndTime;

	void Update()
	{
        if (Input.GetKeyDown(KeyCode.R) && !isRecording)
            StartRecording();

        if (Input.GetKeyDown(KeyCode.T) && isRecording)
            StopRecording();

        if (!isRecording)
            return;

        if (!webCamTexture.didUpdateThisFrame)
            return;

        frameCount++;

        Debug.Log(frameCount);

        frames.Add(webCamTexture.GetPixels32());
	}

    void Start()
    {
        frames = new List<Color32[]>();

        if (webCamResolution == WebCamResolution.Low)
            webCamTexture = new WebCamTexture(WebCamTexture.devices[webCamIndex].name, 160, 120);
        else if (webCamResolution == WebCamResolution.Medium)
            webCamTexture = new WebCamTexture(WebCamTexture.devices[webCamIndex].name, 320, 240);
        else if (webCamResolution == WebCamResolution.High)
            webCamTexture = new WebCamTexture(WebCamTexture.devices[webCamIndex].name, 640, 480);
        else
            webCamTexture = new WebCamTexture(WebCamTexture.devices[webCamIndex].name, 320, 240);

        webCamTexture.Play();

        GetComponent<WebCamTextureRenderer>().GetComponent<GUITexture>().texture = webCamTexture;

        width = webCamTexture.width;
        height = webCamTexture.height;

        bufferSize = width * height * 3;

        buffer = new byte[bufferSize];
    }

    public void StartRecording()
    {
        if (Directory.GetFiles(folderPath, "*", SearchOption.AllDirectories).Length > 0)
            Debug.Log("Target folder is not empty. Recording will not be started.");

        recordingStartTime = Time.realtimeSinceStartup;

        isRecording = true;
    }

    public void StopRecording()
    {
        if (!isRecording)
            return;

        recordingEndTime = Time.realtimeSinceStartup;

        float duration = recordingEndTime - recordingStartTime;

        float fps = frameCount / duration;

        isRecording = false;

        float startTime = Time.realtimeSinceStartup;

        for (int i = 0; i < frames.Count; i++)
        {
            for (int j = 0; j < bufferSize; j += 3)
            {
                buffer[j    ] = frames[i][j / 3].r;
                buffer[j + 1] = frames[i][j / 3].g;
                buffer[j + 2] = frames[i][j / 3].b;
            }

            File.WriteAllBytes(folderPath + "/Frame_" + i, buffer);
        }

        string[] info = new string[3];

        info[0] = "" + 640;
        info[1] = "" + 480;
        info[2] = "" + fps;

        File.WriteAllLines(folderPath + "/info", info);

        Debug.Log(File.ReadAllLines(folderPath + "/info").Length + " lines are read"); 

        float endTime = Time.realtimeSinceStartup;

        Debug.Log((endTime - startTime));
    }
}