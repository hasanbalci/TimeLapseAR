using UnityEngine;
using System.Collections;
using System.Collections.Generic;

[RequireComponent(typeof(ArCamera))]
public class GizmoDrawer : MonoBehaviour
{
    [Header("Gizmos")]
    public bool drawWeightLabels = false;

    public bool drawMarkerEdges = false;
    public Color markerEdgeColor = Color.blue;

    public bool drawControlPoints = false;
    public Color controlPointColor = Color.red;

    public bool drawCorrespondences = false;
    public Color correspondenceColor = Color.green;
    public float pointSize = 3;

    [HideInInspector]
    public Material gizmoMaterial;

    ArCamera arCamera;
    //public ImageIntensityVisualizer imageIntensitVisualizer;

    void Awake()
    {
        arCamera = GameObject.FindObjectOfType<ArCamera>();
    }

    void Start()
    {
        
    }

    void Update()
    {
        if (Input.GetKeyDown(KeyCode.C))
        {
            Application.CaptureScreenshot("C:\\Users\\hasanbalci\\Documents\\MasterThesis\\" + Time.frameCount + ".png");
        }

        if (Input.GetKeyDown(KeyCode.H))
        {
            drawCorrespondences = !drawCorrespondences;
            drawMarkerEdges = !drawMarkerEdges;
            drawControlPoints = !drawControlPoints;
        }
        if (Input.GetKeyDown(KeyCode.L))
        {
            drawCorrespondences = false;
            drawMarkerEdges = false;
            drawControlPoints = false;
        }
    }

    void OnDrawGizmos()
    {
        if(arCamera)
            Gizmos.DrawWireCube(arCamera.bounds.center, arCamera.bounds.size);
    }

    void OnPostRender()
    {
        DrawGizmos();
        //DrawGraph();
    }

    //void DrawGraph()
    //{
    //    lineMaterial.SetPass(0);

    //    GL.LoadPixelMatrix();

    //    int x = 400;
    //    int y = 100;

    //    int yStep = 10;
    //    int xStep = 2;


    //    GL.Begin(GL.LINES);

    //    for (int i = 0; i < arCamera.iterationCounts.Count - 1; i++)
    //    {
    //        Vector2 start = new Vector2(x + i * xStep, y + arCamera.iterationCounts[i] * yStep);
    //        Vector2 end = new Vector2(x + (i + 1) * xStep, y + arCamera.iterationCounts[i + 1] * yStep);
    //        GL.Vertex(start);
    //        GL.Vertex(end);
    //    }

    //    GL.End();
    //}

    //void DrawLine()
    //{
    //    lineMaterial.SetPass(0);

    //    GL.LoadPixelMatrix();

    //    GL.Begin(GL.LINES);

    //    GL.Vertex(imageIntensitVisualizer.start);
    //    GL.Vertex(imageIntensitVisualizer.end);

    //    GL.End();
    //}

    void DrawMatches()
    {
        GL.LoadPixelMatrix();

        List<Vector2> projectedPoints = arCamera.projectedControlPoints;
        List<Vector2> imagePoints = arCamera.imagePoints;
        double[] weights = arCamera.W;

        for (int i = 0; i < projectedPoints.Count; i++)
        {
            //float t = Mathf.InverseLerp((float)arCamera.errMin, (float)arCamera.errMax, (float)arCamera.errcurve[i]);

            gizmoMaterial.color = Color.Lerp(Color.white, correspondenceColor, (float)weights[i]);
            gizmoMaterial.SetPass(0);

            GL.Begin(GL.LINES);
            {
                GL.Vertex(ImageToScreenSpace(projectedPoints[i]));
                GL.Vertex(ImageToScreenSpace(imagePoints[i]));
            }
            GL.End();
        }
    }

    Vector2 ImageToScreenSpace(Vector2 v)
    {
        v.x *= (float)Screen.width / arCamera.width;
        v.y *= (float)Screen.height / arCamera.height;

        return v;
    }

    void DrawMarkerEdges()
    {
        gizmoMaterial.color = markerEdgeColor;
        gizmoMaterial.SetPass(0);

        GL.LoadPixelMatrix();

        GL.Begin(GL.LINES);

        MarkerEdge[] markerEdges = arCamera.markerEdges;

        Camera mainCamera = arCamera.GetComponent<Camera>();

        for (int i = 0; i < markerEdges.Length; i++)
        {
            if (!markerEdges[i].isVisible)
                continue;

            Vector2 s = mainCamera.WorldToScreenPoint(markerEdges[i].visibleStart);
            Vector2 e = mainCamera.WorldToScreenPoint(markerEdges[i].visibleEnd);

            GL.Vertex(s);
            GL.Vertex(e);
        }

        GL.End();
    }

    void OnGUI()
    {
        if (drawWeightLabels)
        {
            List<Vector2> projectedPoints = arCamera.projectedControlPoints;
            double[] weights = arCamera.W;

            for (int i = 0; i < projectedPoints.Count; i++)
            {
                Vector2 p = ImageToScreenSpace(projectedPoints[i]);
                GUI.Label(new Rect(p.x, Screen.height - p.y, 20, 20), "" + weights[i].ToString("0.00"));
            }
        }
    }

    void DrawSamplePoints()
    {
        gizmoMaterial.color = controlPointColor;
        gizmoMaterial.SetPass(0);

        GL.LoadPixelMatrix();

        List<Vector2> projectedPoints = arCamera.projectedControlPoints;
        double[] weights = arCamera.W;


        for (int i = 0; i < projectedPoints.Count; i++)
        {
            //gizmoMaterial.color = Color.Lerp(Color.white, correspondenceColor, (float)weights[i]);
            gizmoMaterial.SetPass(0);


            Vector2 p = projectedPoints[i];

            Vector2 a = p;
            Vector2 b = p;
            Vector2 c = p;
            Vector2 d = p;

            d = ImageToScreenSpace(d);
            c = ImageToScreenSpace(c);
            b = ImageToScreenSpace(b);
            a = ImageToScreenSpace(a);

            float halfPointSize = pointSize / 2;

            a.x -= halfPointSize;
            a.y -= halfPointSize;

            c.x += halfPointSize;
            c.y += halfPointSize;

            d.x -= halfPointSize;
            d.y += halfPointSize;

            b.x += halfPointSize;
            b.y -= halfPointSize;
            GL.Begin(GL.QUADS);
            GL.Vertex(d);
            GL.Vertex(c);
            GL.Vertex(b);
            GL.Vertex(a);
            GL.End();
        }

        
    }

    void DrawGizmos()
    {
        //DrawLine();

        if(drawMarkerEdges)
            DrawMarkerEdges();

        if(drawControlPoints)
            DrawSamplePoints();

        if(drawCorrespondences)
            DrawMatches();
    }
}
