using UnityEngine;
using System.Collections;

public class ControlPoint
{
    public Vector3 position;
    public Vector2 projectedPosition;
    public Vector2 screenNormal;
    public float weight;
    public float diffuseSaliency;
    public int index;
    public int reliability;
    public int dynamicIndex;
    public MarkerEdge markerEdge;
}

public class MarkerEdge
{
    public Vector3 start;
    public Vector3 end;
    public bool isVisible;
    public Vector3 visibleStart;
    public Vector3 visibleEnd;
    public Vector3 n1, n2;
    public bool isHard;
    public bool ignore;
    public bool useIllumination;
    public float intensity;
    public float saliency;
    public ControlPoint[] controlPoints;
};

[DisallowMultipleComponent]
public class Marker : MonoBehaviour
{
    public enum MarkerType
    {
        CadModel,
        ReconstructedMesh
    }

    public MarkerType markerType;

	void Start() 
	{

	}
	
	void Update() 
	{
	
	}
}