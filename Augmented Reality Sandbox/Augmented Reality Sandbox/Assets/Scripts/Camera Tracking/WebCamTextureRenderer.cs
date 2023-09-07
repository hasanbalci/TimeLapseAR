using UnityEngine;
using System.Collections;

public class WebCamTextureRenderer : MonoBehaviour
{
    ArCamera arCamera;

    Rect inset;

    void Awake()
    {
        arCamera = GameObject.FindObjectOfType<ArCamera>();
    }

    void Start()
    {
        inset = new Rect(Screen.width / 2, Screen.height / 2, 0, 0);

        //if (arCamera.usePrerecordedVideo)
        //    guiTexture.texture = arCamera.prerecordedWebcamTexture;
        //else
        //    guiTexture.texture = arCamera.webCamTexture;
    }

    void Update()
    {
        inset = new Rect(Screen.width / 2, Screen.height / 2, 0, 0);
        GetComponent<GUITexture>().pixelInset = inset;
        //guiTexture.texture = tracker.webCamImage;
    }
}
