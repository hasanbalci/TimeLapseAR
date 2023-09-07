using UnityEngine;
//using UnityEditor;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Threading;

public class Cell
{
    public byte[] reliabilityScores;
    public BitArray visibilityMask;
    public Vector3 center;
    public int entropy;
}

[DisallowMultipleComponent]
public class ArCamera : MonoBehaviour
{
    public enum WebCamResolution
    {
        Low,
        Medium,
        High
    }

    #region Public Interface

    public Bounds bounds;

    [Header("Vision Based Tracking")]
    public WebCamResolution webCamResolution = WebCamResolution.High;
    public bool bilinearFiltering = true;
    public float edgeSearchStepSize = 1;
    public float discontinuityThreshold = 3.5f;
    public int edgeSearchRange = 50;
    public int maxMeasurementCount = 1024;
    public bool isTrackingEnabled = false;
    public bool predictCameraMotion = true;
    public bool smoothBorders = true;
    public bool calculateReliabilityInfo = false;
    public bool useReliabilityInfo = false;
    public bool useOnlineReliability = false;
    public bool useVisibilityInfo = false;
    public bool useIlluminationInfo = false;

    [Range(1, 50)]
    [Tooltip("Control point projection interval in screen space (pixels)")]
    public float minProjectedSampleInterval = 20;

    [Header("Marker Processing")]
    public float angleThreshold = 30;
    public float minEdgeLength  = 0;

    [Range(0f, 1f)]
    [Tooltip("Control point generation interval in world space (meters)")]
    public float minSamplingInterval = 0.01f;
    public bool benchmark = false;

    #endregion

    [HideInInspector]
    public bool usePrerecordedVideo;
    [HideInInspector]
    public string folderPath;

    Vector3 translation, rotation;
    Plane[] frustumPlanes;

    [HideInInspector]
    public int webCamIndex;

    [HideInInspector]
    public int width, height;
    Vector2 webCamTextureSize;

    [HideInInspector]
    public WebCamTexture webCamTexture;

    WebCamTextureRenderer webCamTextureRenderer;

    [HideInInspector]
    public List<Vector2> projectedControlPoints;

    [HideInInspector]
    public List<Vector2> imagePoints;

    Color32[] pixels;

    int markerLayerIndex;
    int markerLayerMask;

    [HideInInspector]
    public MarkerEdge[] markerEdges;

    public GameObject combinedMarker;
    Camera[] childCameras;

    int controlPointCount = 0;

    public int totalControlPointCount;
    Camera mainCamera;

    #region Linear Least Squares

    double[,] A;
    double[] B, X;

    [HideInInspector]
    public double[] W;

    #endregion

    #region 3D Grid

    Cell[,,] cells;

    int sizeX, sizeY, sizeZ;

    public float cellSize = 0.025f;

    #endregion

    #region Video Playback

    float recordingFps;

    int frameCount = -1;

    Color32[][] buffers;

    int totalFrameCount;

    #endregion

    #region Benchmarking
    List<Vector3> positions = new List<Vector3>();
    List<Vector3> rotations = new List<Vector3>();
    #endregion

    void Awake()
	{
        SetBounds();
        Init();
        CreateChildCameras();
        CombineMarkers();
        ProcessMarkers();
        InitCells();
        CalculateVisibility();

        if(calculateReliabilityInfo)
            CalculateReliability();

        //ProcessIlluminationInformation();
        //CreateIndicators();
        InitWebCamTexture();
        InitializeSpotLights();
        
	}

    private void SetBounds()
    {
        GameObject g = GameObject.Find("Bounds");
        bounds = g.GetComponent<Collider>().bounds;
        Destroy(g);
    }

    float trackingStartTime;

    void Start()
    {
        testPosition = transform.position;
        testRotation = transform.rotation.eulerAngles;

        if(test)
            SetTestParameters();

        trackingStartTime = Time.realtimeSinceStartup;
    }

    
    Thread[] threads;
    int threadCount;
    void LoadVideo()
    {
        threadCount = SystemInfo.processorCount;
        threads = new Thread[threadCount];
        totalFrameCount = Directory.GetFiles(folderPath, "*", SearchOption.AllDirectories).Length - 1;
        buffers = new Color32[totalFrameCount][];
        float startTime = Time.realtimeSinceStartup;
        for (int i = 0; i < threadCount; i++)
        {
            int threadIndex = i;
            threads[threadIndex] = new Thread(() => LoadVideo_Worker(threadIndex));
            threads[threadIndex].Start();
        }

        for (int i = 0; i < threadCount; i++)
        {
            threads[i].Join();
        }

        float endTime = Time.realtimeSinceStartup;
        
        Debug.Log("Video loaded in: " + (endTime - startTime) + " seconds.");
        Debug.Log("Total frame count: " + totalFrameCount);
    }

    void LoadVideo_Worker(int index)
    {
        int bufferSize = width * height * 3;

        for (int i = index; i < totalFrameCount; i+=threadCount)
        {
            byte[] buffer = File.ReadAllBytes(folderPath + "/Frame_" + i);
            buffers[i] = new Color32[width * height];
            
            for (int j = 0; j < bufferSize; j += 3)
            {
                buffers[i][j / 3].r = buffer[j];
                buffers[i][j / 3].g = buffer[j + 1];
                buffers[i][j / 3].b = buffer[j + 2];
                buffers[i][j / 3].a = 255;
            }
        }
    }

    public Texture2D prerecordedWebcamTexture;

    // init pose tests
    int solveCount = 0;
    int maxSolveCount = 3;
    bool[,] results = new bool[500, 3];
    int type = 0;
    int trialCount = 0;

    [HideInInspector]
    public Vector3 testPosition;
    [HideInInspector]
    public Vector3 testRotation;

    public bool test;

    public float testDistanceThreshold;
    public float testAngleThreshold;

    bool IsPoseCorrect()
    {
        if (Vector3.Distance(transform.position, testPosition) > testDistanceThreshold)
            return false;

        Vector3 rotation = transform.rotation.eulerAngles;

        float angleDistance = Mathf.Pow(Mathf.DeltaAngle(rotation.x, testRotation.x), 2) + Mathf.Pow(Mathf.DeltaAngle(rotation.y, testRotation.y), 2) + Mathf.Pow(Mathf.DeltaAngle(rotation.z, testRotation.z), 2);
        angleDistance = Mathf.Sqrt(angleDistance);

        if (angleDistance > testAngleThreshold)
            return false;

        return true;
    }

    Vector3 ranPos = Vector3.zero;
    Vector3 ranRot = Vector3.zero;

    void SetTestParameters()
    {
        ranPos = testPosition + Random.insideUnitSphere * 0.05f;
        ranRot = testRotation + new Vector3(Random.Range(0f, 1f), Random.Range(0f, 1f), Random.Range(0f, 1f)) * 5F;
    }

    List<float> frameTimes = new List<float>();

    bool pixelProcess = false;

    void Update()
    {
        if (test)
        {
            if (trialCount >= results.GetLength(0))
            {
                //Debug.Log("Done");

                int a = 0;
                int b = 0;
                int c = 0;

                for (int i = 0; i < results.GetLength(0); i++)
                {
                    if (results[i, 0])
                        a++;

                    if (results[i, 1])
                        b++;

                    if (results[i, 2])
                        c++;
                }

                Debug.Log(a + " / " + b + " / " + c);

                return;
            }

            if (solveCount > maxSolveCount)
            {
                results[trialCount, type] = IsPoseCorrect();

                if (type < 2)
                    type++;
                else
                {
                    SetTestParameters();
                    trialCount++;
                    type = 0;
                }

                if (type == 0)
                {
                    useReliabilityInfo = false;
                    useVisibilityInfo = false;
                }
                else if (type == 1)
                {
                    useReliabilityInfo = true;
                    useOnlineReliability = false;
                    useVisibilityInfo = true;
                }
                else if (type == 2)
                {
                    useReliabilityInfo = true;
                    useOnlineReliability = true;
                    useVisibilityInfo = false;
                }

                solveCount = 0;
                transform.position = ranPos;
                transform.rotation = Quaternion.Euler(ranRot);

                return;
            }
        }

        bool shouldUpdateCameraParameters = false;

        if (usePrerecordedVideo)
        {
            if ((int)(Time.realtimeSinceStartup * recordingFps) != frameCount)
            {
                int nextFrameCount = (int)(Time.realtimeSinceStartup * recordingFps);

                if (nextFrameCount - frameCount > 0)
                {
                    frameCount++;

                    int index = (int)Mathf.PingPong(frameCount, totalFrameCount);
                    if (index >= totalFrameCount) index--;

                    pixels = buffers[(int)Mathf.PingPong(index, totalFrameCount)];

                    prerecordedWebcamTexture.SetPixels32(pixels);
                    prerecordedWebcamTexture.Apply();

                    //if (Time.realtimeSinceStartup > 3)
                    //    shouldUpdateCameraParameters = false;
                    //else
                    shouldUpdateCameraParameters = true;
                }

            }

            if (Input.GetKeyDown("p"))
            {
                StartCoroutine(FindPixelFeatures());
            }

            if (Input.GetKeyDown("d"))
            {
                FindLightDirections();
            }
   
            if (Input.GetKeyDown("l"))
            {
                StartCoroutine(FindPixelFeatures());
            }

            if (pixelProcess == true)
            {
                FindLightDirections();
                OptimizeLightValues();
                pixelProcess = false;
            }
            
        }
        else
        {
            if (!webCamTexture.didUpdateThisFrame)
                return;

            webCamTexture.GetPixels32(pixels);

            shouldUpdateCameraParameters = true;
        }

        ProjectPoints();

        if (!isTrackingEnabled)
            return;

        if (shouldUpdateCameraParameters)
        {
            float startTime = Time.realtimeSinceStartup;

            UpdateCameraParameters();
            
            ProjectPoints();
            UpdateCameraParameters();

            float endTime = Time.realtimeSinceStartup;

            ProjectPoints();
            UpdateCameraParameters();

            //ProjectPoints();
            //UpdateCameraParameters();

            //ProjectPoints();
            //UpdateCameraParameters();

            //ProjectPoints();
            //UpdateCameraParameters();

            frameTimes.Add(endTime - startTime);

            positions.Add(transform.position);
            rotations.Add(transform.rotation.eulerAngles);

        }

        //Debug.Log(IsPoseCorrect());
    }

    void OnApplicationQuit()
    {
        float averageFrameTime = 0;

        for (int i = 0; i < totalFrameCount; i++)
        {
            averageFrameTime += frameTimes[i];
        }

        averageFrameTime /= totalFrameCount;

        Debug.Log("Average frame time: " + averageFrameTime);
        
        if (!benchmark)
            return;
        
        System.IO.StreamWriter file = new System.IO.StreamWriter("C:\\Users\\hasanbalci\\Documents\\MasterThesis\\position.txt");

        for (int i = 0; i < positions.Count; i++)
        {
            file.Write(positions[i].x + ";");
            file.Write(positions[i].y + ";");
            file.Write(positions[i].z + "\n");
            //file.Write(rotations[i].x + ";");
            //file.Write(rotations[i].y + ";");
            //file.Write(rotations[i].z + "\n");
        }

        file.Close();

    }

    void Divide(MarkerEdge markerEdge, Vector3 visibleStart, Vector3 visibleEnd, int index, List<ControlPoint> controlPoints)
    {
        if (index >= markerEdge.controlPoints.Length)
            return;

        //use sqrMag here.
        Vector2 edgeDirection = mainCamera.WorldToViewportPoint(visibleEnd) - mainCamera.WorldToViewportPoint(visibleStart);
        edgeDirection.Scale(webCamTextureSize);
        float edgeLength = edgeDirection.magnitude;

        if(edgeLength < minProjectedSampleInterval * 2)
            return;

        ControlPoint controlPoint = markerEdge.controlPoints[index];

        float t = Vector3.Distance(controlPoint.position, markerEdge.start);
        float tMin = Vector3.Distance(visibleStart, markerEdge.start);
        float tMax = Vector3.Distance(visibleEnd, markerEdge.start);

        if(t < tMin)
        {
            Divide(markerEdge, visibleStart, visibleEnd, index * 2 + 2, controlPoints);
        }
        else if(t > tMax)
        {
            Divide(markerEdge, visibleStart, visibleEnd, index * 2 + 1, controlPoints);
        }
        else
        {
            //controlPoint.projectedPosition = mainCamera.WorldToViewportPoint(controlPoint.position);
            //controlPoint.projectedPosition.Scale(webCamTextureSize);
            controlPoints.Add(controlPoint);

            Divide(markerEdge, visibleStart, controlPoint.position, index * 2 + 1, controlPoints);
            Divide(markerEdge, controlPoint.position, visibleEnd, index * 2 + 2, controlPoints);
        }
    }

    List<ControlPoint> visibleControlPoints = new List<ControlPoint>();

    void ProjectPoints()
    {
        frustumPlanes = GeometryUtility.CalculateFrustumPlanes(GetComponent<Camera>());

        float distanceThreshold = 40;
        float angleThreshold = Mathf.Cos(60 * Mathf.Deg2Rad);

        controlPointCount = 0;
        Vector3 cameraPosition = transform.position;

        projectedControlPoints.Clear();
        imagePoints.Clear();
        visibleControlPoints.Clear();

        if (predictCameraMotion)
        {
            transform.Translate(translation * 0.1f);
            transform.Rotate(rotation * 0.1f);
        }

        bool isPrecomputedDataAvailable = bounds.Contains(cameraPosition);

        if (!isPrecomputedDataAvailable) Debug.Log("No data");

        BitArray visibilityMask = null;
        byte[] reliabilityScores = null;


        if(isPrecomputedDataAvailable)
        {
            Vector3 cellIndex = (cameraPosition - bounds.min) / cellSize;
            Cell cell = cells[(int)cellIndex.x, (int)cellIndex.y, (int)cellIndex.z];
            visibilityMask = cell.visibilityMask;
            reliabilityScores = cell.reliabilityScores;
        }

        foreach (MarkerEdge markerEdge in markerEdges)
        {
            ClipEdgeToFrustum(markerEdge);

            if (!markerEdge.isVisible)
                continue;

            Vector3 toEdge = ((markerEdge.end + markerEdge.start) / 2 - cameraPosition).normalized;

            if (Vector3.Dot(toEdge, markerEdge.n1) > 0 && Vector3.Dot(toEdge, markerEdge.n2) > 0)
            {
                //Debug.Log("b");
                continue;
            }

            bool isSilhouette = Vector3.Dot(toEdge, markerEdge.n1) * Vector3.Dot(toEdge, markerEdge.n2) < 0;

            // some raycasts can be avoided if backfacing edges are detected.

            Vector2 edgeDirection = (mainCamera.WorldToViewportPoint(markerEdge.end) - mainCamera.WorldToViewportPoint(markerEdge.start));
            edgeDirection.Scale(webCamTextureSize);
            float edgeLength = edgeDirection.magnitude; // sqr?
            edgeDirection.Normalize();
            Vector2 edgeNormal = new Vector2(-edgeDirection.y, edgeDirection.x);

            List<ControlPoint> controlPoints = new List<ControlPoint>();

            Divide(markerEdge, markerEdge.visibleStart, markerEdge.visibleEnd, 0, controlPoints);

            //int n = Mathf.NextPowerOfTwo((int)(edgeLength / minProjectedSampleInterval)) / 2 - 1;
            //n = Mathf.Min(n, controlPoints.Length);

            for (int i = 0; i < controlPoints.Count; i++)
            {
                if (!TestFrustum(controlPoints[i].position))
                    continue;

                if (useVisibilityInfo && isPrecomputedDataAvailable)
                {
                    if(!visibilityMask[controlPoints[i].index])
                        continue;
                }
                else
                {
                    if (Physics.Linecast(cameraPosition, controlPoints[i].position - (controlPoints[i].position - cameraPosition).normalized * 0.001f, markerLayerMask))
                        continue;
                }

                

                //if (isPrecomputedDataAvailable && useReliabilityInfo)
                //    if (reliabilityScores[controlPoints[i].index] > 4)
                //        continue;

                Vector2 projectedControlPoint = mainCamera.WorldToViewportPoint(controlPoints[i].position);
                projectedControlPoint.Scale(webCamTextureSize);

                if (useReliabilityInfo && useOnlineReliability)
                {
                    controlPoints[i].projectedPosition = projectedControlPoint;
                    controlPoints[i].screenNormal = edgeNormal;
                    controlPoints[i].reliability = 0;

                    MarkerEdge latest = null;

                    for (int j = 0; j < visibleControlPoints.Count; j++)
                    {
                        if (latest == visibleControlPoints[j].markerEdge)
                            continue;

                        if (Vector2.Distance(controlPoints[i].projectedPosition, visibleControlPoints[j].projectedPosition) < distanceThreshold)
                        {
                            if (Mathf.Abs(Vector2.Dot(controlPoints[i].screenNormal, visibleControlPoints[j].screenNormal)) > angleThreshold)
                            {
                                if (controlPoints[i].markerEdge != visibleControlPoints[j].markerEdge)
                                {
                                    controlPoints[i].reliability++;
                                    visibleControlPoints[j].reliability++;
                                    latest = visibleControlPoints[j].markerEdge;
                                }
                            }
                        }
                    }

                    visibleControlPoints.Add(controlPoints[i]);
                    controlPoints[i].dynamicIndex = -1;
                }

                Vector2 targetPoint;
                float intensity;

                if (SearchNearestIntensityDiscontinuity(projectedControlPoint, edgeNormal, discontinuityThreshold, out targetPoint, out intensity))
                {
                    projectedControlPoints.Add(projectedControlPoint);
                    imagePoints.Add(targetPoint);

                    float distance = Vector3.Dot((targetPoint - projectedControlPoint), edgeNormal);

                    B[controlPointCount] = distance;

                    for (int k = 0; k < childCameras.Length; k++)
                    {
                        Vector2 point = childCameras[k].WorldToViewportPoint(controlPoints[i].position);
                        point.Scale(webCamTextureSize);

                        A[controlPointCount, k] = Vector2.Dot(point - projectedControlPoint, edgeNormal);
                    }

                    float w = 1;

                    //if(!isSilhouette && useIlluminationInfo)
                    //    w = controlPoints[i].diffuseSaliency;

                    if (!isSilhouette)
                        w *= 0.75f;


                    Bounds safeBounds = bounds;
                    safeBounds.min += new Vector3(cellSize, cellSize, cellSize) * 0.5f;
                    safeBounds.max -= new Vector3(cellSize, cellSize, cellSize) * 0.5f;

                    if (isPrecomputedDataAvailable && useReliabilityInfo && safeBounds.Contains(cameraPosition) && !useOnlineReliability)
                    {
                        w = w / (GetReliabilityScoreTrilinear(controlPoints[i], cameraPosition) + 1);
                    }

                    //w = (edgeSearchRange * 4 - distance);

                    if (smoothBorders)
                    {
                        float minDistX = Mathf.Min(projectedControlPoint.x - 0, width - projectedControlPoint.x);
                        float minDistY = Mathf.Min(projectedControlPoint.y - 0, height - projectedControlPoint.y);
                        float minDistToBounds = Mathf.Min(minDistX, minDistY);

                        if(minDistToBounds < 25) w *= (minDistToBounds / 25);
                    }

                    controlPoints[i].dynamicIndex = controlPointCount;

                    W[controlPointCount++] = w;
                }
            }
        }

        if (useReliabilityInfo && useOnlineReliability)
        {
            for (int i = 0; i < visibleControlPoints.Count; i++)
            {
                if (visibleControlPoints[i].dynamicIndex != -1)
                    W[visibleControlPoints[i].dynamicIndex] /= (visibleControlPoints[i].reliability + 1);
            }
        }

        //Debug.Log(controlPointCount);
    }

    bool TestFrustum(Vector3 point)
    {
        for (int i = 0; i < frustumPlanes.Length; i++)
        {
            if (frustumPlanes[i].GetDistanceToPoint(point) <= 0)
                return false;
        }

        return true;
    }

    bool ClipEdgeToFrustum(MarkerEdge edge)
    {
        bool isStartPointInside = TestFrustum(edge.start);
        bool isEndPointInside = TestFrustum(edge.end);

        edge.isVisible = true;

        float length = Vector3.Distance(edge.start, edge.end);

        if (isStartPointInside && isEndPointInside)
        {
            edge.isVisible = true;
            edge.visibleStart = edge.start;
            edge.visibleEnd = edge.end;
            return true;
        }
        else if (isStartPointInside && !isEndPointInside)
        {
            Ray ray = new Ray(edge.start, edge.end - edge.start);
            float enter;
            float minDistance = float.MaxValue;

            for (int i = 0; i < frustumPlanes.Length; i++)
            {
                frustumPlanes[i].Raycast(ray, out enter);

                if (enter > 0 && enter < minDistance)
                {
                    minDistance = enter;
                }
            }

            edge.visibleStart = edge.start;
            edge.visibleEnd = edge.start + ray.direction * minDistance;

            return true;
        }
        else if (!isStartPointInside && isEndPointInside)
        {
            Ray ray = new Ray(edge.end, edge.start - edge.end);
            float enter;
            float minDistance = float.MaxValue;

            for (int i = 0; i < frustumPlanes.Length; i++)
            {
                frustumPlanes[i].Raycast(ray, out enter);

                if (enter > 0 && enter < minDistance)
                {
                    minDistance = enter;
                }
            }

            edge.visibleEnd = edge.end;
            edge.visibleStart = edge.end + ray.direction * minDistance;

            return true;
        }
        else
        {
            List<float> distances = new List<float>(8);

            Ray ray = new Ray(edge.start, edge.end - edge.start);
            float enter;

            for (int i = 0; i < frustumPlanes.Length; i++)
            {
                frustumPlanes[i].Raycast(ray, out enter);

                if (enter > 0 && enter < length)
                {
                    distances.Add(enter);
                }
            }

            distances.Sort();

            for (int i = 1; i < distances.Count; i++)
            {
                float t = (distances[i - 1] + distances[i]) / 2;

                Vector3 point = ray.origin + t * ray.direction;

                if (TestFrustum(point))
                {
                    edge.visibleStart = ray.origin + distances[i - 1] * ray.direction;
                    edge.visibleEnd = ray.origin + distances[i] * ray.direction;
                    return true;
                }
            }
        }

        edge.isVisible = false;

        return false;
    }

    void UpdateCameraParameters()
    {
        solveCount++;

        int info;
        alglib.lsfitreport report;

        float startTime = Time.realtimeSinceStartup;

        try
        {
            alglib.lsfitlinearw(B, W, A, controlPointCount, 6, out info, out X, out report);
        }
        catch (System.Exception e)
        {
            return;
        }
        
        float endTime = Time.realtimeSinceStartup;

        //Debug.Log("Solved in " + (endTime - startTime) + " seconds. (n = " + controlPointCount + ")");

        //Debug.Log(X[0] + " " + X[1] + " " + X[2] + " " + X[3] + " " + X[4] + " " + X[5]);

        translation = new Vector3((float)X[0], (float)X[1], (float)X[2]) * 0.1f;
        rotation = new Vector3((float)X[3], (float)X[4], (float)X[5]);

        transform.Translate(translation);
        transform.Rotate(rotation);
    }

    bool SearchNearestIntensityDiscontinuity(Vector2 origin, Vector2 direction, float intensityThreshold, out Vector2 result, out float intensity)
    {
        Vector2 offset = Vector2.zero;
        direction *= edgeSearchStepSize;
        intensity = 0;

        origin -= (direction / 2);

        Vector2 position = origin;

        if (!(position.x > 1 && position.x < width - 1 && position.y > 1 && position.y < height - 1))
        {
            result = Vector2.zero;
            return false;
        }

        float vPos = GetImageIntensity(position);
        float vNeg = vPos;
        float v;

        float lastNegDifference = 0;
        float lastPosDifference = 0;

        for (int i = 0; i < edgeSearchRange; i++)
        {
            offset += direction;
            position = origin + offset;

            if (position.x > 1 && position.x < width - 1 && position.y > 1 && position.y < height - 1)
            {
                v = GetImageIntensity(position);
                float intensityDifference = Mathf.Abs(v - vPos);

                if (lastPosDifference > intensityThreshold && intensityDifference < lastPosDifference)
                {
                    result = position - direction * 1.5f;

                    intensity = intensityDifference;

                    return true;
                }

                lastPosDifference = intensityDifference;
                vPos = v;
            }

            position = origin - offset;

            if (position.x > 1 && position.x < width - 1 && position.y > 1 && position.y < height - 1)
            {
                v = GetImageIntensity(position);
                float intensityDifference = Mathf.Abs(v - vNeg);

                if (lastNegDifference > intensityThreshold && intensityDifference < lastNegDifference)
                {
                    result = position + direction * 1.5f;

                    intensity = intensityDifference;

                    return true;
                }

                lastNegDifference = intensityDifference;
                vNeg = v;
            }
        }

        result = Vector2.zero;

        return false;
    }

    public float GetImageIntensity(Vector2 position)
    {
        if (bilinearFiltering)
        {
            position.x = position.x - 0.5f;
            position.y = position.y - 0.5f;

            int x = (int)position.x;
            int y = (int)position.y;

            float u_ratio = position.x - x;
            float v_ratio = position.y - y;

            float u_opposite = 1 - u_ratio;

            int b = x + y * width;

            return (ToGrayscale(pixels[b]) * u_opposite + ToGrayscale(pixels[b + 1]) * u_ratio) * (1 - v_ratio) +
                   (ToGrayscale(pixels[b + width]) * u_opposite + ToGrayscale(pixels[b + width + 1]) * u_ratio) * v_ratio;
        }
        else
        {
            return ToGrayscale(pixels[(int)position.x + (int)position.y * width]);
        }
    }

    int ToGrayscale(Color32 color)
    {
        return (color.r + color.g + color.b) / 3;
    }

    void Init()
    {
        A = new double[maxMeasurementCount, 6];
        X = new double[6];
        B = new double[maxMeasurementCount];
        W = new double[maxMeasurementCount];

        projectedControlPoints = new List<Vector2>();
        imagePoints = new List<Vector2>();

        webCamTextureRenderer = GameObject.FindObjectOfType<WebCamTextureRenderer>();
        mainCamera = GetComponent<Camera>();
    }

    void InitWebCamTexture()
    {
        if (usePrerecordedVideo)
        {
            string[] info = File.ReadAllLines(folderPath + "/info");

            width = int.Parse(info[0]);
            height = int.Parse(info[1]);
            recordingFps = float.Parse(info[2]);
        }
        else
        {
            if (webCamResolution == WebCamResolution.Low)
                webCamTexture = new WebCamTexture(WebCamTexture.devices[webCamIndex].name, 160, 120);
            else if (webCamResolution == WebCamResolution.Medium)
                webCamTexture = new WebCamTexture(WebCamTexture.devices[webCamIndex].name, 320, 240);
            else if (webCamResolution == WebCamResolution.High)
                webCamTexture = new WebCamTexture(WebCamTexture.devices[webCamIndex].name, 640, 480);
            else
                webCamTexture = new WebCamTexture(WebCamTexture.devices[webCamIndex].name, 320, 240);

            webCamTexture.Play();

            width = webCamTexture.width;
            height = webCamTexture.height;
        }

        webCamTextureSize.x = width;
        webCamTextureSize.y = height;

        Debug.Log(width + " " + height);

        pixels = new Color32[width * height];

        prerecordedWebcamTexture = new Texture2D(width, height);

        if (usePrerecordedVideo)
            webCamTextureRenderer.GetComponent<GUITexture>().texture = prerecordedWebcamTexture;
        else
            webCamTextureRenderer.GetComponent<GUITexture>().texture = webCamTexture;

        if (usePrerecordedVideo)
            LoadVideo();
    }

    void CreateChildCameras()
    {
        childCameras = new Camera[6];

        for (int i = 0; i < 3; i++)
        {
            GameObject tCamera = new GameObject("T" + ((char)('X' + i)).ToString());
            tCamera.AddComponent<Camera>().fieldOfView = GetComponent<Camera>().fieldOfView;
            tCamera.transform.parent = transform;
            Vector3 position = Vector3.zero;
            position[i] = 0.1f;
            tCamera.transform.localPosition = position;
            tCamera.transform.localRotation = Quaternion.Euler(0, 0, 0);
            tCamera.SetActive(false);
            //tCamera.hideFlags |= HideFlags.HideInHierarchy;
            childCameras[i] = tCamera.GetComponent<Camera>();

            GameObject rCamera = new GameObject("R" + ((char)('X' + i)).ToString());
            rCamera.AddComponent<Camera>().fieldOfView = GetComponent<Camera>().fieldOfView;
            rCamera.transform.parent = transform;
            Vector3 rotation = Vector3.zero;
            rotation[i] = 1;
            rCamera.transform.localPosition = Vector3.zero;
            rCamera.transform.localRotation = Quaternion.Euler(rotation);
            rCamera.SetActive(false);
            childCameras[i + 3] = rCamera.GetComponent<Camera>();
            //rCamera.hideFlags |= HideFlags.HideInHierarchy;
        }
    }

    public void CombineMarkers()
    {
        if (combinedMarker)
            DestroyImmediate(combinedMarker);

        int freeLayerIndex = -1;

        for (int i = 0; i < 32; i++)
        {
            if (LayerMask.LayerToName(i).Equals(""))
            {
                freeLayerIndex = i;
                break;
            }
        }

        if (freeLayerIndex == -1)
        {
            Debug.Log("No free layer is found.");
        }

        markerLayerIndex = freeLayerIndex;

        markerLayerMask = 1 << markerLayerIndex;

        //Debug.Log("Layer " + freeLayerIndex + " is reserved.");

        Marker[] markers = GameObject.FindObjectsOfType<Marker>();

        if (markers.Length == 0)
        {
            Debug.Log("No marker is found.");
            return;
        }

        //List<MeshFilter> meshFilters = new List<MeshFilter>();

        //foreach (Marker marker in markers)
        //{
        //    meshFilters.AddRange(marker.GetComponentsInChildren<MeshFilter>());
        //}

        //CombineInstance[] combine = new CombineInstance[meshFilters.Count];

        //for (int i = 0; i < meshFilters.Count; i++)
        //{
        //    combine[i].mesh = meshFilters[i].sharedMesh;
        //    combine[i].transform = meshFilters[i].transform.localToWorldMatrix;
        //    meshFilters[i].gameObject.SetActive(false);
        //}

        //combinedMarker = new GameObject("Combined Marker");



        //MeshFilter meshFilter = combinedMarker.AddComponent<MeshFilter>();
        //meshFilter.sharedMesh = new Mesh();
        //meshFilter.sharedMesh.CombineMeshes(combine, false);
        //combinedMarker.AddComponent<MeshCollider>();

        

        combinedMarker = Instantiate(markers[0].gameObject) as GameObject;
        combinedMarker.name = "Combined Marker";
        combinedMarker.AddComponent<MeshCollider>();
        combinedMarker.GetComponent<Renderer>().enabled = false;
        combinedMarker.layer = markerLayerIndex;
        for (int i = 0; i < 32; i++)
            Physics.IgnoreLayerCollision(markerLayerIndex, i, true);

        markers[0].gameObject.SetActive(false);
    }

    public void ProcessMarkers()
    {
        totalControlPointCount = 0;

        float overlapThreshold = float.Epsilon;
        List<MarkerEdge> extractedEdges = new List<MarkerEdge>();

        MeshFilter mf = combinedMarker.GetComponentInChildren<MeshFilter>();

        int triangleCount = mf.sharedMesh.triangles.Length / 3;
        int[] triangles = mf.sharedMesh.triangles;
        Vector3[] vertices = mf.sharedMesh.vertices;

        for (int i = 0; i < vertices.Length; i++)
            vertices[i] = mf.transform.TransformPoint(vertices[i]);

        List<Vector3> normals = new List<Vector3>();
        List<MarkerEdge> edges = new List<MarkerEdge>();

        for (int i = 0; i < triangleCount; ++i)
        {
            Vector3 v1, v2;
            v1 = vertices[triangles[3 * i + 1]] - vertices[triangles[3 * i]];
            v2 = vertices[triangles[3 * i + 2]] - vertices[triangles[3 * i]];

            normals.Add(Vector3.Cross(v1, v2).normalized);
            
        }

        int submeshCount = mf.sharedMesh.subMeshCount;
        //Debug.Log("Submesh count:" + submeshCount);

        int[] triangleIDs = new int[triangleCount];

        int offset = 0;

        for (int i = 0; i < submeshCount; i++)
        {
            int submeshLength = mf.sharedMesh.GetIndices(i).Length / 3;

            for (int j = 0; j < submeshLength; j++)
                triangleIDs[j + offset] = i;

            offset += submeshLength; 
        }

        // Loop of shame. Should optimize this sometime.
        for (int i = 0; i < triangleCount; ++i)
        {
            for (int j = i + 1; j < triangleCount; ++j)
            {
                for (int x = 0; x < 3; ++x)
                {
                    for (int y = 0; y < 3; ++y)
                    {
                        int s1, e1, s2, e2;

                        s1 = triangles[i * 3 + x];
                        e1 = triangles[i * 3 + (x + 1) % 3];

                        s2 = triangles[j * 3 + y];
                        e2 = triangles[j * 3 + (y + 1) % 3];

                        bool overlap = false;

                        if ((vertices[s1] - vertices[e2]).magnitude < overlapThreshold && (vertices[s2] - vertices[e1]).magnitude < overlapThreshold)
                            overlap = true;

                        if (overlap)
                        {
                            float dot = Vector3.Dot(normals[i], normals[j]);

                            if (Mathf.Acos(dot) > angleThreshold * Mathf.Deg2Rad && Vector3.Distance(vertices[s1], vertices[e1]) > minEdgeLength)
                            {
                                MarkerEdge edge = new MarkerEdge();
                                edge.start = vertices[s1];
                                edge.end = vertices[e1];

                                edge.n1 = normals[i];
                                edge.n2 = normals[j];

                                edge.useIllumination = (triangleIDs[i] == triangleIDs[j]);

                                edge.isHard = true;

                                GenerateControlPoints(edge);

                                edges.Add(edge);
                            }
                        }
                    }
                }
            }
        }

        extractedEdges.AddRange(edges);

        markerEdges = extractedEdges.ToArray();

        Debug.Log(markerEdges.Length + " edges are extracted from the CAD model.");
        Debug.Log(totalControlPointCount + " control points are generated.");
    }

    void GenerateControlPoints(MarkerEdge markerEdge)
    {
        float edgeLength = Vector3.Distance(markerEdge.start, markerEdge.end);

        int edgeCount = 1;
        int controlPointCount = 1;
        int subdivisionDepth = 0;

        Vector3 edgeDirection = (markerEdge.end - markerEdge.start).normalized;

        List<ControlPoint> controlPointList = new List<ControlPoint>();

        while (edgeLength > 2 * minSamplingInterval)
        {
            float controlPointOffset = edgeLength / 2;
            float controlPointStep = edgeLength;

            for (int i = 0; i < edgeCount; i++)
            {
                ControlPoint controlPoint = new ControlPoint();
                controlPoint.index = totalControlPointCount++;
                controlPoint.position = markerEdge.start + (controlPointOffset + i * controlPointStep) * edgeDirection;
                controlPoint.markerEdge = markerEdge;
                controlPointList.Add(controlPoint);
            }

            subdivisionDepth++;
            controlPointCount += edgeCount;
            edgeLength /= 2;
            edgeCount *= 2;
        }

        markerEdge.controlPoints = controlPointList.ToArray();
    }

    void ProcessIlluminationInformation()
    {
        Light[] lights = GameObject.FindObjectsOfType<Light>();

        foreach (MarkerEdge markerEdge in markerEdges)
        {
            
            foreach (ControlPoint controlPoint in markerEdge.controlPoints)
            {
                if (markerEdge.useIllumination)
                {
                    controlPoint.diffuseSaliency = 0;

                    foreach (Light light in lights)
                    {
                        Vector3 toLightNormalized = (light.transform.position - controlPoint.position).normalized;
                        controlPoint.diffuseSaliency += Mathf.Abs(Mathf.Max(0, Vector3.Dot(markerEdge.n1, toLightNormalized)) - Mathf.Max(0, Vector3.Dot(markerEdge.n2, toLightNormalized)));
                    }
                }
                else
                    controlPoint.diffuseSaliency = 1;
            }
        }
    }

    public GameObject markerEdgeModel, controlPointModel;

    void CreateIndicators()
    {
        foreach (MarkerEdge markerEdge in markerEdges)
        {
            GameObject g = GameObject.Instantiate(markerEdgeModel, (markerEdge.start + markerEdge.end) / 2, Quaternion.identity) as GameObject;

            float edgeLength = Vector3.Distance(markerEdge.start, markerEdge.end);

            g.transform.LookAt(markerEdge.end);
            g.transform.localScale = new Vector3(0.001f, 0.001f, edgeLength);
            g.GetComponent<Renderer>().material.color = Color.blue;

            foreach (ControlPoint controlPoint in markerEdge.controlPoints)
            {
                GameObject g2 = GameObject.Instantiate(controlPointModel, controlPoint.position, Quaternion.identity) as GameObject;
                g2.transform.localScale = new Vector3(0.002f, 0.002f, 0.002f);
                g2.GetComponent<Renderer>().material.color = Color.Lerp(Color.black, Color.white, controlPoint.diffuseSaliency);
            }
        }
    }

    public void InitCells()
    {
        sizeX = (int)(bounds.size.x / cellSize) + 1;
        sizeY = (int)(bounds.size.y / cellSize) + 1;
        sizeZ = (int)(bounds.size.z / cellSize) + 1;

        cells = new Cell[sizeX, sizeY, sizeZ];

        Debug.Log("3D grid dimensions: " + sizeX + " " + sizeY + " " + sizeZ);

        for (int i = 0; i < sizeX; i++)
        {
            for (int j = 0; j < sizeY; j++)
            {
                for (int k = 0; k < sizeZ; k++)
                {
                    Cell cell = new Cell();
                    cell.visibilityMask = new BitArray(totalControlPointCount);
                    cell.reliabilityScores = new byte[totalControlPointCount];
                    cell.center = new Vector3((i + 0.5f) * cellSize, (j + 0.5f) * cellSize, (k + 0.5f) * cellSize) + bounds.min;
                    cells[i, j, k] = cell;
                }
            }
        }
    }

    public void CalculateVisibility()
    {
        float startTime = Time.realtimeSinceStartup;

        foreach(Cell cell in cells)
        {
            foreach (MarkerEdge markerEdge in markerEdges)
            {
                foreach (ControlPoint controlPoint in markerEdge.controlPoints)
                {
                    Vector3 start = cell.center;
                    Vector3 end = controlPoint.position - (controlPoint.position - start).normalized * 0.001f;

                    cell.visibilityMask[controlPoint.index] = !Physics.Linecast(start, end, markerLayerMask); 
                }
            }
        }

        //for (int i = 0; i < 1000; i++)
        //{
        //    Cell cell = cells[Random.Range(0, cells.GetLength(0)), Random.Range(0, cells.GetLength(1)), Random.Range(0, cells.GetLength(2))];
        //    Cell otherCell = cells[Random.Range(0, cells.GetLength(0)), Random.Range(0, cells.GetLength(1)), Random.Range(0, cells.GetLength(2))];

        //    BitArray result = cell.visibilityMask.Xor(otherCell.visibilityMask);

        //    bool equals = true;

        //    for(int j = 0; j < result.Count; j++)
        //    {
        //        if (result[j])
        //        {
        //            equals = false;
        //            break;
        //        }
        //    }

        //    if (equals)
        //        Debug.Log("Equals");
        //}

        float endTime = Time.realtimeSinceStartup;

        Debug.Log("Visibility information is calculated in " + (endTime - startTime).ToString("0.000") + " seconds.");
    }

    public void CalculateReliability()
    {
        threadCount = SystemInfo.processorCount;
        threads = new Thread[threadCount];

        float startTime = Time.realtimeSinceStartup;

        for (int i = 0; i < threadCount; i++)
        {
            int threadIndex = i;
            threads[threadIndex] = new Thread(() => CalculateReliability_Worker(threadIndex));
            threads[threadIndex].Start();
        }

        for (int i = 0; i < threadCount; i++)
        {
            threads[i].Join();
        }


        float endTime = Time.realtimeSinceStartup;

        Debug.Log("Reliability information is calculated in " + (endTime - startTime).ToString("0.000") + " seconds.");
    }

    public void CalculateReliability_Worker(int index)
    {
        float angleThreshold = Mathf.Cos(2f * Mathf.Deg2Rad);
        float dotThreshold = Mathf.Cos(60 * Mathf.Deg2Rad);

        for(int i = 0; i < cells.GetLength(0); i++)
        {
            for(int j = 0; j < cells.GetLength(1); j++)
            {
                for (int k = index; k < cells.GetLength(2); k += threadCount)
                {
                    Cell cell = cells[i, j, k];

                    foreach (MarkerEdge markerEdge in markerEdges)
                    {
                        Vector3 edgeDirection = (markerEdge.end - markerEdge.start);

                        foreach (ControlPoint controlPoint in markerEdge.controlPoints)
                        {
                            if (!cell.visibilityMask[controlPoint.index])
                                continue;

                            Vector3 toControlPoint = (controlPoint.position - cell.center).normalized;
                            Vector3 normal = Vector3.Cross(toControlPoint, edgeDirection).normalized;

                            foreach (MarkerEdge otherMarkerEdge in markerEdges)
                            {
                                if (markerEdge == otherMarkerEdge)
                                    continue;

                                Vector3 otherEdgeDirection = (otherMarkerEdge.end - otherMarkerEdge.start);

                                float prevAngle = float.NegativeInfinity;

                                foreach (ControlPoint otherControlPoint in otherMarkerEdge.controlPoints)
                                {
                                    if (!cell.visibilityMask[otherControlPoint.index])
                                        continue;

                                    Vector3 toOtherControlPoint = (otherControlPoint.position - cell.center).normalized;


                                    float angle = Vector3.Dot(toControlPoint, toOtherControlPoint);

                                    if (angle > prevAngle)
                                    {
                                        prevAngle = angle;
                                    }
                                    else break;

                                    if (angle > angleThreshold)
                                    {
                                        Vector3 otherNormal = Vector3.Cross(toOtherControlPoint, otherEdgeDirection).normalized;

                                        if (Mathf.Abs(Vector3.Dot(normal, otherNormal)) > dotThreshold)
                                        {
                                            cell.reliabilityScores[controlPoint.index]++;
                                            
                                        }

                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }

    float GetReliabilityScoreTrilinear(ControlPoint controlPoint, Vector3 cameraPosition)
    {
        int index = controlPoint.index;

        Vector3 position = ((cameraPosition - bounds.min) / cellSize) - Vector3.one * 0.5f;

        int x = (int)position.x;
        int y = (int)position.y;
        int z = (int)position.z;

        float u_ratio = position.x - x;
        float v_ratio = position.y - y;
        float w_ratio = position.z - z;

        float u_opposite = 1 - u_ratio;
        float v_opposite = 1 - v_ratio;
        float w_opposite = 1 - w_ratio;

        float score = (((cells[x, y, z].reliabilityScores[index] * u_opposite + cells[x + 1, y, z].reliabilityScores[index] * u_ratio) * v_opposite +
                      (cells[x, y + 1, z].reliabilityScores[index] * u_opposite + cells[x + 1, y + 1, z].reliabilityScores[index] * u_ratio) * v_ratio)) * w_opposite +

                      (((cells[x, y, z + 1].reliabilityScores[index] * u_opposite) + cells[x + 1, y, z + 1].reliabilityScores[index] * u_ratio) * v_opposite +
                      (cells[x, y + 1, z + 1].reliabilityScores[index] * u_opposite + cells[x + 1, y + 1, z + 1].reliabilityScores[index] * u_ratio) * v_ratio) * w_ratio;

        return score;
    }

    void OnGUI()
    {
        //GUI.Label(new Rect(10, 10, 100, 20), "" + frameCount);
    }

    void OnDrawGizmos()
    {
        //if (bounds != null)
            //Gizmos.DrawWireCube(bounds.center, bounds.size);
    }

    public Vector3 center;
    public float radius;
    public List<GameObject> spotlight = new List<GameObject>(8);
    public GameObject pointlight;

    void InitializeSpotLights()
    {
        //center = new Vector3(-0.575f, 0.0f, 0.09f);
        //radius = 0.5f;
        spotlight[0].transform.position = center + new Vector3(0, radius * (Mathf.Sin(Mathf.PI / 6)), radius * (Mathf.Cos(Mathf.PI / 6)));
        spotlight[1].transform.position = center + new Vector3(0, radius * (Mathf.Sin(Mathf.PI / 3)), radius * (Mathf.Cos(Mathf.PI / 3)));
        spotlight[2].transform.position = center + new Vector3(0, radius * (Mathf.Sin(Mathf.PI / 3)), -radius * (Mathf.Cos(Mathf.PI / 3)));
        spotlight[3].transform.position = center + new Vector3(0, radius * (Mathf.Sin(Mathf.PI / 6)), -radius * (Mathf.Cos(Mathf.PI / 6)));
        spotlight[4].transform.position = center + new Vector3(-radius * (Mathf.Cos(Mathf.PI / 6)), radius * (Mathf.Sin(Mathf.PI / 6)), 0);
        spotlight[5].transform.position = center + new Vector3(-radius * (Mathf.Cos(Mathf.PI / 3)), radius * (Mathf.Sin(Mathf.PI / 3)), 0);
        spotlight[6].transform.position = center + new Vector3(radius * (Mathf.Cos(Mathf.PI / 3)), radius * (Mathf.Sin(Mathf.PI / 3)), 0);
        spotlight[7].transform.position = center + new Vector3(radius * (Mathf.Cos(Mathf.PI / 6)), radius * (Mathf.Sin(Mathf.PI / 6)), 0);

        for (int i = 0; i < spotlight.Count; i++ )
            spotlight[i].transform.rotation = Quaternion.LookRotation(center - spotlight[i].transform.position);

    }

    List<Vector3> pixelPoints = new List<Vector3>();
    List<Vector3> pixelNormals = new List<Vector3>();

    byte[] rgbImageBuffer;

    Texture2D rgbImage;

    Color32[] rgbPixelsAll;
    List<int> shadingPixelsAll = new List<int>();
    List<int> shadingPixelsSamples = new List<int>();

    IEnumerator FindPixelFeatures()
    {    
        Ray ray;
        RaycastHit hit;

        Application.CaptureScreenshot("C:/Users/hasanbalci/Documents/MasterThesis/TestImage.png");
        yield return new WaitForFixedUpdate();

        rgbImageBuffer = File.ReadAllBytes("C:/Users/hasanbalci/Documents/MasterThesis/TestImage.png");        
        rgbImage = new Texture2D(width, height);

        rgbImage.LoadImage(rgbImageBuffer);
        rgbPixelsAll = rgbImage.GetPixels32();

        for (int i = 0; i < rgbPixelsAll.Length; i++)
            shadingPixelsAll.Add(ToGrayscale(rgbPixelsAll[i]));

        for (int i = 0; i < height; i = i + 1)
        {
            for (int j = 0; j < width; j = j + 1)
            {
                ray = Camera.main.ScreenPointToRay(new Vector3(j, i, 0));
                if (Physics.Raycast(ray, out hit))
                {
                    pixelPoints.Add(hit.point);
                    pixelNormals.Add(hit.normal);
                    shadingPixelsSamples.Add(shadingPixelsAll[i * width + j]);
                }
            }
        }

        pixelProcess = true;

        //Debug.Log(pixelPoints.Count);
        Debug.Log("Pixel points are found.");
        Debug.Log("Pixel normals are found.");
        Debug.Log("Shading pixel values are found.");
    }

    List<List<int>> allLightDirectionsBoolean = new List<List<int>>();
    List<List<Vector3>> allLightDirectionsValue = new List<List<Vector3>>();

    void FindLightDirections()
    {
        Vector3 direction;
        RaycastHit hit;
        List<int> lightDirectionsPerPixelBoolean = new List<int>();
        List<Vector3> lightDirectionsPerPixelValue = new List<Vector3>();

        for (int i = 0; i < pixelPoints.Count; i++)
        {
            for (int j = 0; j < spotlight.Count; j++)
            {
                direction = pixelPoints[i] - spotlight[j].transform.position;
                if (Physics.Raycast(spotlight[j].transform.position, direction, out hit))
                {
                    if (hit.point == pixelPoints[i])
                    {
                        lightDirectionsPerPixelBoolean.Add(1);
                        lightDirectionsPerPixelValue.Add(direction);
                    }
                    else
                    {
                        lightDirectionsPerPixelBoolean.Add(0);
                        lightDirectionsPerPixelValue.Add(new Vector3(0, 0, 0));
                    }
                }
                else
                {
                    lightDirectionsPerPixelBoolean.Add(0);
                    lightDirectionsPerPixelValue.Add(new Vector3(0, 0, 0));
                }
            }

            allLightDirectionsBoolean.Add(lightDirectionsPerPixelBoolean);
            allLightDirectionsValue.Add(lightDirectionsPerPixelValue);

            lightDirectionsPerPixelBoolean = new List<int>();
            lightDirectionsPerPixelValue = new List<Vector3>();
        }

        Debug.Log("Light directions are found.");      

    }

    void FunctionVector(double[] lightValues, double[] functions, object obj)
    {
        double estimatedPixelValue = 0;

        for (int i = 0; i < pixelPoints.Count; i++)
        {
            for (int j = 0; j < spotlight.Count; j++)
            {
                if (allLightDirectionsBoolean[i][j] == 1)
                    estimatedPixelValue += Vector3.Dot((float)lightValues[j] * allLightDirectionsValue[i][j], pixelNormals[i]);
            }

            estimatedPixelValue += lightValues[spotlight.Count];

            functions[i] = shadingPixelsSamples[i] - estimatedPixelValue;
        }
    }

    public GameObject ball1;
    public GameObject ball2;
    public GameObject ball3;
    public GameObject cube;

    void OptimizeLightValues()
    {
        double[] lightValues = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
        double epsg = 0.0000000001;
        double epsf = 0;
        double epsx = 0;
        int maxits = 0;
        alglib.minlmstate state;
        alglib.minlmreport report;

        alglib.minlmcreatev(pixelPoints.Count, lightValues, 0.0001, out state);
        alglib.minlmsetcond(state, epsg, epsf, epsx, maxits);
        alglib.minlmoptimize(state, FunctionVector, null, null);
        alglib.minlmresults(state, out lightValues, out report);
        
        Debug.Log(alglib.ap.format(lightValues, 4));
        Debug.Log(report.terminationtype);
        Debug.Log(report.iterationscount);
        Debug.Log("Light values are found.");

        for (int i = 0; i < spotlight.Count; i++)
            spotlight[i].GetComponent<Light>().intensity = (float)lightValues[i] * 20;

        pointlight.GetComponent<Light>().intensity = (float)lightValues[8] * 20;
        //pointlight.GetComponent<Light>().intensity = 0.1f;

        ball1.transform.position = new Vector3(-0.456f, 0.020f, 0.089f);
        cube.transform.position = new Vector3(-0.64f, 0.020f, 0.106f);
        ball3.transform.position = new Vector3(-0.404f, 0.020f, 0.167f);

    }

}